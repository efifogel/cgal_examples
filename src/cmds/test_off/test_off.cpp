#include <iostream>
#include <string>
#include <vector>
#include <list>

#define USE_SURFACE_MESH

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#if defined(USE_SURFACE_MESH)
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#else
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/draw_polyhedron.h>
#endif

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>

#if ! defined(USE_SURFACE_MESH)
#include "Extended_polyhedron_items.h"
#endif

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[]) {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

#if defined(USE_SURFACE_MESH)
  using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
  using Face_color_map = Mesh::Property_map<Mesh::Face_index, CGAL::IO::Color>;
#else
  using Traits = CGAL::Polyhedron_traits_with_normals_3<Kernel>;
  using Mesh = CGAL::Polyhedron_3<Traits, Extended_polyhedron_items>;
  using Face_color_map = boost::property_map<Mesh, CGAL::dynamic_face_property_t<CGAL::IO::Color> >::type;
#endif

  using Vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
  using Edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
  using Face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

  const char* filename = (argc > 1) ? argv[1] : "cube.off";
  std::cout << filename << std::endl;
  Mesh mesh;

  bool has_colors = false;
  Face_color_map fcolors;

#if defined(USE_SURFACE_MESH)
  if (! CGAL::IO::read_polygon_mesh(filename, mesh, CGAL::parameters::verbose(true).repair_polygon_soup(true))) {
    std::cerr << "Error: could not read mesh from '" << filename << "'.\n";
    return 1;
  }

  auto face_property_maps = mesh.properties<Face_descriptor>();
  for (const auto& pm : face_property_maps) std::cout << pm << std::endl;
#else
  if (! CGAL::IO::read_polygon_mesh(filename, mesh, CGAL::parameters::verbose(true).repair_polygon_soup(true).
                                    face_color_map(fcolors))) {
    std::cerr << "Error: could not read mesh from '" << filename << "'.\n";
    return 1;
  }
#endif
  std::cout << "Loaded mesh with " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces.\n";

  CGAL::Graphics_scene_options<Mesh, Vertex_descriptor, Edge_descriptor, Face_descriptor> gso;
  gso.ignore_all_vertices(true);
  gso.ignore_all_edges(true);
  gso.colored_face = [](const Mesh&, typename boost::graph_traits<Mesh>::face_descriptor) -> bool
  { return true; };
  gso.face_color =  [] (const Mesh&, typename boost::graph_traits<Mesh>::face_descriptor fh) -> CGAL::IO::Color {
    if (fh == boost::graph_traits<Mesh>::null_face()) return CGAL::IO::Color(100, 125, 200);
    return get_random_color(CGAL::get_default_random());
  };

#if defined(USE_SURFACE_MESH)
  auto fcolors_optional = mesh.property_map<Face_descriptor, CGAL::IO::Color>("f:color");
  if (fcolors_optional) {
    has_colors = true;
    fcolors = fcolors_optional.value();
  }
  else std::cerr << "Warning: could not read colors from '" << filename << "'.\n";
#else
  CGAL::IO::Color zero_color(0, 0, 0, 255);
  for (auto fd : CGAL::faces(mesh)) {
    CGAL::IO::Color c = get(fcolors, fd);
    if (c != zero_color) has_colors = true;
  }
#endif

  if (has_colors) {
    std::size_t idx = 0;
    for (auto fd : CGAL::faces(mesh)) {
      CGAL::IO::Color c = get(fcolors, fd);
      std::cout << "face[" << idx++ << "]: ("
                << int(c.red())   << ", "
                << int(c.green()) << ", "
                << int(c.blue())  << ", "
                << int(c.alpha()) << ")\n";
    }
    gso.colored_face=[](const Mesh&, Face_descriptor) -> bool { return true; };
    gso.face_color=[&](const Mesh& mesh, Face_descriptor fd) -> CGAL::IO::Color { return get(fcolors, fd); };
  }

  CGAL::draw(mesh, gso, filename);

  // General validity
  auto is_valid = mesh.is_valid();
  if (! is_valid) std::cerr << "The mesh is not valid\n";

  // Closed
  auto is_closed = CGAL::is_closed(mesh);
  if (! is_closed) std::cerr << "The mesh is not closed\n";

  // Triangular mesh
  auto is_tri = CGAL::is_triangle_mesh(mesh);
  if (! is_tri) std::cerr << "The mesh is not triangular\n";

  // Does self intersect
  auto self_intersect = PMP::does_self_intersect(mesh);
  if (self_intersect) std::cerr << "The mesh self intersects\n";

  // Connected components
#if defined(USE_SURFACE_MESH)
  auto fccmap = mesh.add_property_map<Face_descriptor, std::size_t>("f:CC").first;
  std::size_t num_ccs = PMP::connected_components(mesh, fccmap);
  std::cout << "Number of connected components: " << num_ccs << "\n";

  // Does bound a volume
  if (is_closed && is_tri) {
    std::vector<bool> isoo(num_ccs);
    bool bound_volume =
      PMP::does_bound_a_volume(mesh, CGAL::parameters::is_cc_outward_oriented(std::reference_wrapper(isoo)));
    if (! bound_volume) std::cerr << "The mesh does not bound a volume\n";
    for (auto ccid = 0; ccid < num_ccs; ++ccid)
      if (! isoo[ccid]) std::cerr << "Component " << ccid << " is not outward oriented\n";
  }
#endif

  // if (num_ccs > 1) {
  //   std::vector<Mesh> ccs;
  //   CGAL::Polygon_mesh_processing::split_connected_components(mesh, ccs);
  //   assert(ccs.size() == 2);
  //   std::size_t i = 0;
  //   for (auto& cc : ccs) {
  //     CGAL::draw(cc, gso, filename);
  //     bool does_bound_volume = PMP::does_bound_a_volume(cc);
  //     std::cout << "The mesh [" << i++ << "]" << ((does_bound_volume) ? "does" : "does not") << " bound a volume\n";
  //   }
  // }
  return 0;
}
