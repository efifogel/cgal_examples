#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include "Extended_polyhedron_items.h"

int main(int argc, char* argv[]) {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Traits = CGAL::Polyhedron_traits_with_normals_3<Kernel>;
  using Mesh = CGAL::Polyhedron_3<Traits, Extended_polyhedron_items>;
  Kernel kernel;
  const char* filename = (argc > 1) ? argv[1] : "cube.off";
  std::cout << filename << std::endl;
  Mesh mesh;
  auto rc = CGAL::IO::read_polygon_mesh(filename, mesh, CGAL::parameters::verbose(true));
  CGAL::draw(mesh, filename);
  return 0;
}
