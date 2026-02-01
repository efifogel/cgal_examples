#include <CGAL/Bbox_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include "cgalex/merge_coplanar_faces.h"
#include "cgalex/corefine_difference.h"

#include "Extended_polyhedron_items.h"
#include "Face_normal_map.h"
#include "gaussian_maps.h"
#include "minkowski_sum_3.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = Kernel::Point_3;
using Traits = CGAL::Polyhedron_traits_with_normals_3<Kernel>;
using Mesh = CGAL::Polyhedron_3<Traits, Extended_polyhedron_items>;

/*! Convert a bounding box to a polhedron representation.
 */
Mesh bbox_to_polyhedron(const CGAL::Bbox_3& bbox) {
  std::vector<Point_3> points(8);
  points[0] = Point_3(bbox.xmin(), bbox.ymin(), bbox.zmin());
  points[1] = Point_3(bbox.xmax(), bbox.ymin(), bbox.zmin());
  points[2] = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmin());
  points[3] = Point_3(bbox.xmin(), bbox.ymax(), bbox.zmin());
  points[4] = Point_3(bbox.xmin(), bbox.ymin(), bbox.zmax());
  points[5] = Point_3(bbox.xmax(), bbox.ymin(), bbox.zmax());
  points[6] = Point_3(bbox.xmax(), bbox.ymax(), bbox.zmax());
  points[7] = Point_3(bbox.xmin(), bbox.ymax(), bbox.zmax());

  Mesh mesh;
  CGAL::convex_hull_3(points.begin(), points.end(), mesh);
  return mesh;
}

/*! Main entry.
 */
int main(int argc, char* argv[]) {
  Kernel kernel;
  const char* filename1 = (argc > 1) ? argv[1] : "star.off";
  const char* filename2 = (argc > 2) ? argv[2] : "cube.off";

  // Construct the first mesh
  Mesh mesh1;
  auto rc1 = CGAL::IO::read_polygon_mesh(filename1, mesh1);
  CGAL::draw(mesh1, filename1);

  // Compute the bounding box
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh1);
  auto bbox_mesh1 = bbox_to_polyhedron(bbox);
  std::cout << bbox << std::endl;

  // Compute the difference between the bounding box and the mesh
  Mesh complement_mesh1;
  // corefine_difference(bbox_mesh1, mesh1, complement_mesh1);

  // Decompose the difference into convex pieces and compute their Gaussian maps
  std::vector<Arrangement> gms1;
  gaussian_maps(mesh1, std::back_inserter(gms1), kernel);
  std::cout << filename1 << " decomposed into " << gms1.size() << " pieces\n";

  Mesh mesh2;
  auto rc2 = CGAL::IO::read_polygon_mesh(filename2, mesh2);
  CGAL::draw(mesh2, filename2);
  std::vector<Arrangement> gms2;
  gaussian_maps(mesh2, std::back_inserter(gms2), kernel);
  std::cout << filename2 << " decomposed into " << gms2.size() << " pieces\n";
  auto ms = minkowski_sum_3<Mesh>(gms1.begin(), gms1.end(), gms2.begin(), gms2.end());
  auto np = CGAL::parameters::geom_traits(kernel);
  Face_normal_map<Mesh> normals;
  CGAL::Polygon_mesh_processing::compute_face_normals(ms, normals, np);
  merge_coplanar_faces(ms, normals, np);
  CGAL::draw(ms, "Minkowski Sum");

  return 0;
}
