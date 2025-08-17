#include <iostream>
#include <string>

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

#if ! defined(USE_SURFACE_MESH)
#include "Extended_polyhedron_items.h"
#endif

int main(int argc, char* argv[]) {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

#if defined(USE_SURFACE_MESH)
  using Mesh   = CGAL::Surface_mesh<Kernel::Point_3>;
#else
  using Traits = CGAL::Polyhedron_traits_with_normals_3<Kernel>;
  using Mesh = CGAL::Polyhedron_3<Traits, Extended_polyhedron_items>;
  using Face_color_map = boost::property_map<Mesh, CGAL::dynamic_face_property_t<CGAL::IO::Color> >::type;
#endif

  using Face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

  const char* filename = (argc > 1) ? argv[1] : "cube.off";
  std::cout << filename << std::endl;
  Mesh mesh;

#if defined(USE_SURFACE_MESH)
  if (! CGAL::IO::read_polygon_mesh(filename, mesh, CGAL::parameters::verbose(true).repair_polygon_soup(true))) {
    std::cerr << "Error: could not read mesh from '" << filename << "'.\n";
    return 1;
  }
  std::cout << "Loaded mesh with " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces.\n";

  auto face_property_maps = mesh.properties<Face_descriptor>();
  for (const auto& pm : face_property_maps) std::cout << pm << std::endl;

  // Show colors for faces
  auto fcolors_optional = mesh.property_map<Face_descriptor, CGAL::IO::Color >("f:color");
  if (! fcolors_optional) {
     std::cerr << "Error: could not read colors from '" << filename << "'.\n";
     return 1;
  }
  auto fcolors = fcolors_optional.value();

#else
  Face_color_map fcolors;
  if (! CGAL::IO::read_polygon_mesh(filename, mesh, CGAL::parameters::verbose(true).repair_polygon_soup(true).
                                    face_color_map(fcolors))) {
    std::cerr << "Error: could not read mesh from '" << filename << "'.\n";
    return 1;
  }
  std::cout << "Loaded mesh with " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces.\n";
#endif

  std::size_t idx = 0;
  for (auto fd : CGAL::faces(mesh)) {
    CGAL::IO::Color c = get(fcolors, fd);
    std::cout << "face[" << idx++ << "]: ("
              << int(c.red())   << ", "
              << int(c.green()) << ", "
              << int(c.blue())  << ", "
              << int(c.alpha()) << ")\n";
  }

  CGAL::draw(mesh, filename);

  return 0;
}
