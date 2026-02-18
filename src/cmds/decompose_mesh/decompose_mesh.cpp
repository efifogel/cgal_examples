#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

namespace PMP = CGAL::Polygon_mesh_processing;

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Vector_3 = typename Kernel::Vector_3;
using Polyhedral_mesh = CGAL::Polyhedron_3<Kernel>;

int main(int argc, char* argv[]) {
  const char* filename = (argc > 1) ? argv[1] : "cd_bug.off";
  Polyhedral_mesh mesh;

  CGAL::IO::read_polygon_mesh(filename, mesh);

  if (mesh.is_empty()) std::cerr << "The mesh is empty\n";
  if (! mesh.is_valid()) std::cerr << "The mesh is invalid\n";
  if (! CGAL::is_closed(mesh)) std::cerr << "The mesh is open\n";
  if (PMP::does_self_intersect(mesh)) std::cerr << "The mesh self intersects\n";
  if (! PMP::is_outward_oriented(mesh)) std::cerr << "The mesh is not outward oriented\n";
  using Ed = typename boost::graph_traits<Polyhedral_mesh>::edge_descriptor;
  std::vector<Ed> edges;
  PMP::degenerate_edges(mesh, std::back_inserter(edges));
  using Fd = typename boost::graph_traits<Polyhedral_mesh>::face_descriptor;
  if (! edges.empty()) std::cerr << "The mesh has degenerate edges\n";
  std::vector<Fd> faces;
  PMP::degenerate_faces(mesh, std::back_inserter(faces));
  if (! faces.empty()) std::cerr << "The mesh has degenerate faces\n";

  CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> mesh_nef(mesh);
  CGAL::convex_decomposition_3(mesh_nef);
  return 0;
}
