#include <iostream>
#include <fstream>
#include <string>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point = Kernel::Point_3;
using Plane = Kernel::Plane_3;
using Vector = Kernel::Vector_3;
using Mesh = CGAL::Surface_mesh<Point>;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[]) {
  const char* filename = (argc > 1) ? argv[1] : "cd_bug.off";
  Mesh mesh;
  if (! CGAL::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Error: Cannot read input.off" << std::endl;
    return -1;
  }

  if (mesh.is_empty()) std::cerr << "The mesh is empty\n";
  if (! mesh.is_valid()) std::cerr << "The mesh is invalid\n";
  if (! CGAL::is_closed(mesh)) std::cerr << "The mesh is open\n";
  if (! PMP::is_outward_oriented(mesh)) std::cerr << "The mesh is not outward oriented\n";
  if (! PMP::is_outward_oriented(mesh)) std::cerr << "The mesh is not outward oriented\n";
  using Ed = typename boost::graph_traits<Mesh>::edge_descriptor;
  std::vector<Ed> edges;
  PMP::degenerate_edges(mesh, std::back_inserter(edges));
  using Fd = typename boost::graph_traits<Mesh>::face_descriptor;
  if (! edges.empty()) std::cerr << "The mesh has degenerate edges\n";
  std::vector<Fd> faces;
  PMP::degenerate_faces(mesh, std::back_inserter(faces));
  if (! faces.empty()) std::cerr << "The mesh has degenerate faces\n";

  // --- CALCULATE MIDPOINT ---
  CGAL::Bbox_3 bbox = PMP::bbox(mesh);
  std::cout << "Choose axis to split (x, y, or z): ";
  char axis;
  std::cin >> axis;

  Point pivot;
  Vector normal;
  double mid;

  if (axis == 'x' || axis == 'X') {
    mid = (bbox.xmin() + bbox.xmax()) / 2.0;
    pivot = Point(mid, 0, 0);
    normal = Vector(1, 0, 0);
  }
  else if (axis == 'y' || axis == 'Y') {
    mid = (bbox.ymin() + bbox.ymax()) / 2.0;
    pivot = Point(0, mid, 0);
    normal = Vector(0, 1, 0);
  }
  else {
    mid = (bbox.zmin() + bbox.zmax()) / 2.0;
    pivot = Point(0, 0, mid);
    normal = Vector(0, 0, 1);
    axis = 'z'; // Default to z if input is weird
  }

  std::cout << "Splitting at " << axis << " = " << mid << "..." << std::endl;

  Plane plane_neg(pivot, normal);
  Plane plane_pos(pivot, -normal);

  // Perform the Split
  Mesh mesh_pos = mesh; // "Positive" side of the plane
  Mesh mesh_neg = mesh; // "Negative" side of the plane

  // Keep the "inner" side (opposite to normal)
  // To get the part BELOW the midpoint, use the UPward normal
  PMP::clip(mesh_neg, plane_pos, PMP::parameters::clip_volume(true));

  // To get the part ABOVE the midpoint, use the DOWNward normal
  PMP::clip(mesh_pos, plane_neg, PMP::parameters::clip_volume(true));

  // Export
  CGAL::IO::write_polygon_mesh("part_positive.off", mesh_pos, CGAL::parameters::stream_precision(17));
  CGAL::IO::write_polygon_mesh("part_negative.off", mesh_neg, CGAL::parameters::stream_precision(17));

  std::cout << "Exported part_positive.off and part_negative.off" << std::endl;

  return 0;
}
