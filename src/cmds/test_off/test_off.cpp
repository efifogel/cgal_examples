#include <iostream>
#include <string>
#include <vector>
#include <list>

#define USE_SURFACE_MESH

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#if defined(USE_SURFACE_MESH)
#include <CGAL/Surface_mesh.h>
#if defined(CGALEX_HAS_VISUAL)
#include <CGAL/draw_surface_mesh.h>
#endif
#else
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#if defined(CGALEX_HAS_VISUAL)
#include <CGAL/draw_polyhedron.h>
#endif
#endif

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#if ! defined(USE_SURFACE_MESH)
#include "Extended_polyhedron_items.h"
#endif

#include "cgalex/retriangulate_faces.h"

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

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

  using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
  using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
  using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
  using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

  const char* filename = (argc > 1) ? argv[1] : "cube.off";
  // std::cout << filename << std::endl;
  Mesh mesh;

  bool has_colors = false;
  Face_color_map fcolors;

#if defined(USE_SURFACE_MESH)
  if (! CGAL::IO::read_polygon_mesh(filename, mesh, params::verbose(true).repair_polygon_soup(true))) {
    std::cerr << "Error: could not read mesh from '" << filename << "'.\n";
    return 1;
  }

  auto face_property_maps = mesh.properties<face_descriptor>();
  // for (const auto& pm : face_property_maps) std::cout << pm << std::endl;
#else
  if (! CGAL::IO::read_polygon_mesh(filename, mesh, params::verbose(true).repair_polygon_soup(true).
                                    face_color_map(fcolors))) {
    std::cerr << "Error: could not read mesh from '" << filename << "'.\n";
    return 1;
  }
#endif
  // std::cout << "Loaded mesh with " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces.\n";

#if defined(CGALEX_HAS_VISUAL)
  CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso;
  gso.ignore_all_vertices(true);

  gso.ignore_all_edges(true);
  gso.colored_edge = [](const Mesh& mesh, typename boost::graph_traits<Mesh>::edge_descriptor ed) -> bool {
    auto hd = CGAL::halfedge(ed, mesh);
    auto ohd = CGAL::opposite(hd, mesh);
    if ((CGAL::face(hd, mesh) == boost::graph_traits<Mesh>::null_face()) ||
        (CGAL::face(ohd, mesh) == boost::graph_traits<Mesh>::null_face())) {
      std::cout << "ed: " << ed << std::endl;
      return true;
    }
    return false;
  };
  gso.edge_color =  [] (const Mesh&, typename boost::graph_traits<Mesh>::edge_descriptor ed) -> CGAL::IO::Color {
    return CGAL::IO::Color(0, 128, 0);
  };

  gso.colored_face = [](const Mesh&, typename boost::graph_traits<Mesh>::face_descriptor) -> bool { return true; };
  gso.face_color =  [] (const Mesh&, typename boost::graph_traits<Mesh>::face_descriptor fh) -> CGAL::IO::Color {
    if (fh == boost::graph_traits<Mesh>::null_face()) return CGAL::IO::Color(0, 0, 0);
    return get_random_color(CGAL::get_default_random());
    // return CGAL::IO::Color(128, 0, 0);
  };

#if defined(USE_SURFACE_MESH)
  auto fcolors_optional = mesh.property_map<face_descriptor, CGAL::IO::Color>("f:color");
  if (fcolors_optional) {
    has_colors = true;
    fcolors = fcolors_optional.value();
  }
  // else std::cerr << "Warning: could not read colors from '" << filename << "'.\n";
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
    gso.colored_face=[](const Mesh&, face_descriptor) -> bool { return true; };
    gso.face_color=[&](const Mesh& mesh, face_descriptor fd) -> CGAL::IO::Color { return get(fcolors, fd); };
  }
#endif

  Kernel kernel;

  // Repair
  // PMP::stitch_borders(mesh, params::apply_per_connected_component(true));
  PMP::remove_degenerate_faces(mesh);
  PMP::remove_degenerate_edges(mesh);
  // CGAL::draw(mesh, gso, filename);

  // General validity
  auto is_valid = mesh.is_valid();
  if (! is_valid) std::cerr << "The mesh is not valid\n";

  // Closed
  auto is_closed = CGAL::is_closed(mesh);
  if (! is_closed) {
    std::cerr << "The mesh is not closed\n";
    std::vector<halfedge_descriptor> border_cycles;
    PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
    std::cout << "# boundary cycles: " << border_cycles.size() << std::endl;
    std::cout << "Cycle: " << std::endl;
    for (const auto& cycle : border_cycles) std::cout << cycle << std::endl;
  }

  // Does self intersect
  auto self_intersect = PMP::does_self_intersect(mesh);
  if (self_intersect) std::cerr << "The mesh self intersects\n";

  // // Detect Zero area
  // std::vector<face_descriptor> zero_area_faces;
  // // const double area_threshold = 1e-15; // Define a small threshold for floating point comparison
  // const double area_threshold = 0;
  // for (auto f : mesh.faces()) {
  //   auto area = CGAL::Polygon_mesh_processing::face_area(f, mesh);
  //   if (area <= area_threshold) zero_area_faces.push_back(f);
  // }
  // if (! zero_area_faces.empty()) {
  //   std::cerr << "The mesh has " << zero_area_faces.size() << " zero-area faces\n";
  //   auto res = PMP::remove_almost_degenerate_faces(mesh, params::geom_traits(kernel));
  //   if (! res) std::cout << "Almost degenerate faces have been removed\n";
  // }

  if (is_closed) {
    using Vector_3 = typename Kernel::Vector_3;
    auto np = params::geom_traits(kernel);
    auto normals = mesh.add_property_map<face_descriptor, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;

#if 1
    PMP::compute_face_normals(mesh, normals, np);
    retriangulate_faces(mesh, normals, np);
#else
    Mesh temp;
    PMP::remesh_planar_patches(mesh, temp);
    std::swap(mesh, temp);
#endif

    // General validity
    is_valid = mesh.is_valid(true);
    if (! is_valid) std::cerr << "The mesh is not valid 2\n";

    // Closed
    is_closed = CGAL::is_closed(mesh);
    if (! is_closed) {
      std::cerr << "The mesh is not closed 2\n";
      std::vector<halfedge_descriptor> border_cycles;
      PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));
      std::cout << "# boundary cycles: " << border_cycles.size() << std::endl;
      std::cout << "Cycle: " << std::endl;
      auto vpm = get_property_map(CGAL::vertex_point, mesh);
      for (const auto& cycle : border_cycles) {
        auto v_src = CGAL::source(cycle, mesh);
        const auto& src = get(vpm, v_src);
        auto v_trg = CGAL::target(cycle, mesh);
        const auto& trg = get(vpm, v_trg);
        std::cout << src << " " << trg << std::endl;
      }
    }

    // Does self intersect
    auto self_intersect = PMP::does_self_intersect(mesh);
    if (self_intersect) std::cerr << "The mesh self intersects\n";

#if defined(CGALEX_HAS_VISUAL)
    CGAL::draw(mesh, gso, "Intermediate");
#endif
  }

  // Triangular mesh
  auto is_tri = CGAL::is_triangle_mesh(mesh);
  if (! is_tri) {
    std::cerr << "The mesh is not triangular; triangulating\n";
    // PMP::triangulate_faces(mesh, params::geom_traits(kernel));
    // PMP::compute_face_normals(mesh, normals, np);
    // triangulate_faces(mesh, normals);
    // CGAL::draw(mesh, gso, "Final");

    // // General validity
    // auto is_valid = mesh.is_valid();
    // if (! is_valid) std::cerr << "The mesh is not valid\n";

    // // Closed
    // auto is_closed = CGAL::is_closed(mesh);
    // if (! is_closed) std::cerr << "The mesh is not closed\n";
  }

  // Connected components
#if defined(USE_SURFACE_MESH)
  auto fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
  std::size_t num_ccs = PMP::connected_components(mesh, fccmap);
  std::cout << "Number of connected components: " << num_ccs << "\n";

  // Does bound a volume
  if (is_closed && is_tri) {
    std::vector<bool> isoo(num_ccs);
    bool bound_volume =
      PMP::does_bound_a_volume(mesh, params::is_cc_outward_oriented(std::reference_wrapper(isoo)));
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
