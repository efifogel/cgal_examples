#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/draw_arrangement_2.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_traits_with_normals_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/convex_decomposition_3.h>

#include "cgalex/merge_coplanar_faces.h"
#include "cgalex/corefine_difference.h"
#include "cgalex/corefine_union.h"
#include "cgalex/retriangulate_faces.h"
#include "cgalex/soup_difference.h"
#include "cgalex/soup_union.h"

#include "Adder.h"
#include "arr_gaussian_map.h"
#include "Extended_polyhedron_items.h"
#include "Face_normal_map.h"
#include "gaussian_map.h"
#include "Polyhedron_builder.h"

/*! Splits a mesh into two halves using a specified plane.
 *
 * \param input_mesh The source mesh to split.
 * \param split_plane The plane used to cut the mesh.
 * \param[out] out_mesh_pos The resulting mesh on the positive side of the normal.
 * \param[out] out_mesh_neg The resulting mesh on the negative side of the normal.
 * \param close_holes If true, caps the cut area to make watertight volumes.
 * \return true if the split was successful.
 */
template <typename Mesh_, typename Kernel_>
bool split_mesh_by_plane(const Mesh_& input_mesh, const typename Kernel_::Plane_3& split_plane,
                         Mesh_& mesh_pos, Mesh_& mesh_neg, bool close_holes = true) {
  // Check if input is valid
  if (CGAL::is_empty(input_mesh)) {
    std::cerr << "Input mesh is empty!" << std::endl;
    return false;
  }

  // Copy the input mesh into both outputs because clip() modifies in-place.
  mesh_pos = input_mesh;
  mesh_neg = input_mesh;

  // Compute the Negative Side
  // PMP::clip() removes the volume on the POSITIVE side of the plane,
  // leaving us with the negative side.
  bool neg_success = PMP::clip(mesh_neg, split_plane, params::clip_volume(close_holes));

  if (! neg_success) {
    std::cerr << "Failed to clip for the negative side." << std::endl;
    return false;
  }

  // Compute the Positive Side
  // To get the positive side, we clip using the OPPOSITE plane.
  // This removes the "negative" side relative to the original plane.
  bool pos_success = PMP::clip(mesh_pos, split_plane.opposite(), params::clip_volume(close_holes));

  if (! pos_success) {
    std::cerr << "Failed to clip for the positive side." << std::endl;
    return false;
  }

  // Garbage collection
  // PMP operations might leave isolated vertices or deleted elements marked.
  mesh_pos.collect_garbage();
  mesh_neg.collect_garbage();

  return true;
}

// // Example usage
// int main() {
//     Mesh mesh;

//     // Create a simple cube for testing (corners from -1 to 1)
//     CGAL::make_hexahedron(
//         Point(-1,-1,-1), Point(1,-1,-1), Point(1,1,-1), Point(-1,1,-1),
//         Point(-1,-1,1),  Point(1,-1,1),  Point(1,1,1),  Point(-1,1,1),
//         mesh
//     );

//     // Define a Vertical Plane (e.g., X = 0)
//     // Plane equation: 1*x + 0*y + 0*z + 0 = 0
//     Plane vertical_plane(1, 0, 0, 0);

//     Mesh positive_part, negative_part;

//     std::cout << "Splitting mesh at X=0..." << std::endl;

//     if (split_mesh_by_plane(mesh, vertical_plane, positive_part, negative_part)) {
//         std::cout << "Success!" << std::endl;
//         std::cout << "Right Half (Pos) Vertices: " << positive_part.number_of_vertices() << std::endl;
//         std::cout << "Left Half (Neg) Vertices:  " << negative_part.number_of_vertices() << std::endl;

//         // Save to file
//         CGAL::IO::write_polygon_mesh("split_pos.off", positive_part);
//         CGAL::IO::write_polygon_mesh("split_neg.off", negative_part);
//     } else {
//         std::cerr << "Split operation failed." << std::endl;
//     }

//     return 0;
// }

template <typename Kernel>
typename Kernel::Iso_cuboid_3 make_iso_cube(double size, const Kernel& kernel) {
  size *= 0.5;
  typename Kernel::Point_3 p(-size, -size, -size);
  typename Kernel::Point_3 q(size, size, size);
  return typename Kernel::Iso_cuboid_3(p, q);
}

/*! Convert a bounding box to a polhedron representation.
 */
template <typename Mesh_, typename Cube_>
Mesh_ to_mesh(const Cube_& bbox, double offset = 0.0) {
  using Mesh = Mesh_;
  using Point_pm = typename boost::property_map<Mesh, boost::vertex_point_t>::type;
  using Point_3 = typename boost::property_traits<Point_pm>::value_type;
  std::vector<Point_3> points(8);
  points[0] = Point_3(bbox.xmin()-offset, bbox.ymin()-offset, bbox.zmin()-offset);
  points[1] = Point_3(bbox.xmax()+offset, bbox.ymin()-offset, bbox.zmin()-offset);
  points[2] = Point_3(bbox.xmax()+offset, bbox.ymax()+offset, bbox.zmin()-offset);
  points[3] = Point_3(bbox.xmin()-offset, bbox.ymax()+offset, bbox.zmin()-offset);
  points[4] = Point_3(bbox.xmin()-offset, bbox.ymin()-offset, bbox.zmax()+offset);
  points[5] = Point_3(bbox.xmax()+offset, bbox.ymin()-offset, bbox.zmax()+offset);
  points[6] = Point_3(bbox.xmax()+offset, bbox.ymax()+offset, bbox.zmax()+offset);
  points[7] = Point_3(bbox.xmin()-offset, bbox.ymax()+offset, bbox.zmax()+offset);

  Mesh mesh;
  CGAL::convex_hull_3(points.begin(), points.end(), mesh);
  return mesh;
}

/*! Convert a bounding box to a polhedron representation.
 */
template <typename Mesh_, typename Cube_>
bool split_mesh(const Cube_& bbox, Mesh_& mesh_pos, Mesh_& mesh_neg, double offset = 0.0) {
  using Mesh = Mesh_;
  using Point_pm = typename boost::property_map<Mesh, boost::vertex_point_t>::type;
  using Point_3 = typename boost::property_traits<Point_pm>::value_type;
  auto midx = (bbox.xmin() + bbox.xmax()) / 2;

  std::vector<Point_3> points_pos(8);
  points_pos[0] = Point_3(bbox.xmin()-offset, bbox.ymin()-offset, bbox.zmin()-offset);
  points_pos[1] = Point_3(midx+offset, bbox.ymin()-offset, bbox.zmin()-offset);
  points_pos[2] = Point_3(midx+offset, bbox.ymax()+offset, bbox.zmin()-offset);
  points_pos[3] = Point_3(bbox.xmin()-offset, bbox.ymax()+offset, bbox.zmin()-offset);
  points_pos[4] = Point_3(bbox.xmin()-offset, bbox.ymin()-offset, bbox.zmax()+offset);
  points_pos[5] = Point_3(midx+offset, bbox.ymin()-offset, bbox.zmax()+offset);
  points_pos[6] = Point_3(midx+offset, bbox.ymax()+offset, bbox.zmax()+offset);
  points_pos[7] = Point_3(bbox.xmin()-offset, bbox.ymax()+offset, bbox.zmax()+offset);
  CGAL::convex_hull_3(points_pos.begin(), points_pos.end(), mesh_pos);

  std::vector<Point_3> points_neg(8);
  points_neg[0] = Point_3(midx-offset, bbox.ymin()-offset, bbox.zmin()-offset);
  points_neg[1] = Point_3(bbox.xmax()+offset, bbox.ymin()-offset, bbox.zmin()-offset);
  points_neg[2] = Point_3(bbox.xmax()+offset, bbox.ymax()+offset, bbox.zmin()-offset);
  points_neg[3] = Point_3(midx-offset, bbox.ymax()+offset, bbox.zmin()-offset);
  points_neg[4] = Point_3(midx-offset, bbox.ymin()-offset, bbox.zmax()+offset);
  points_neg[5] = Point_3(bbox.xmax()+offset, bbox.ymin()-offset, bbox.zmax()+offset);
  points_neg[6] = Point_3(bbox.xmax()+offset, bbox.ymax()+offset, bbox.zmax()+offset);
  points_neg[7] = Point_3(midx-offset, bbox.ymax()+offset, bbox.zmax()+offset);
  CGAL::convex_hull_3(points_neg.begin(), points_neg.end(), mesh_neg);

  return true;
}

/*!
 */
template <typename Mesh, typename OutputIterator, typename Kernel>
OutputIterator gaussian_maps(Mesh& mesh, OutputIterator oi, const Kernel& kernel) {
  if (CGAL::is_strongly_convex_3(mesh, kernel)) {
    *oi++ = gaussian_map<Arrangement>(mesh, kernel);
    return oi;
  }
  CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> mesh_nef(mesh);
  CGAL::convex_decomposition_3(mesh_nef);
  auto ci = ++mesh_nef.volumes_begin();
  for (; ci != mesh_nef.volumes_end(); ++ci) {
    if (! ci->mark()) continue;
    Mesh mesh;
    mesh_nef.convert_inner_shell_to_polyhedron(ci->shells_begin(), mesh);
    draw(mesh, "Convex part");
    *oi++ = gaussian_map<Arrangement>(mesh, kernel);
  }
  return oi;
}

/*!
 */
template <typename Mesh, typename InputIterator, typename Arrangement_>
Mesh minkowski_sum_3(InputIterator begin, InputIterator end, const Arrangement_& gm2) {
  CGAL::Arr_face_overlay_traits<Arrangement, Arrangement, Arrangement, Adder> overlay_traits;
  Mesh ms;
  bool first = true;
  for (auto it = begin; it != end; ++it) {
    Arrangement gm;
    CGAL::overlay(*it, gm2, gm, overlay_traits); /* \label{lst:minkSum:overlay} */
    using Halfedge_ds = typename Mesh::HalfedgeDS;
    Polyhedron_builder<Halfedge_ds, Arrangement> surface(gm); /* \label{lst:minkSum:surface} */
    if (first) {
      ms.delegate(surface);
      PMP::triangulate_faces(ms);
      first = false;
      continue;
    }
    Mesh mesh;
    mesh.delegate(surface);
    PMP::triangulate_faces(mesh);
    soup_union(ms, mesh, ms);
  }
  return ms;
}

/*! Main entry.
 */
int main(int argc, char* argv[]) {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  using Point_3 = Kernel::Point_3;
  using Vector_3 = typename Kernel::Vector_3;
  using Traits = CGAL::Polyhedron_traits_with_normals_3<Kernel>;
  using Polyhedral_mesh = CGAL::Polyhedron_3<Traits, Extended_polyhedron_items>;
  using Surface_mesh = CGAL::Surface_mesh<Point_3>;
  using Fd = typename boost::graph_traits<Surface_mesh>::face_descriptor;

  Kernel kernel;
  auto np = params::geom_traits(kernel);

  const char* filename = (argc > 1) ? argv[1] : "star.off";

  // Construct the first mesh
  Surface_mesh mesh_sm;
  auto rc1 = CGAL::IO::read_polygon_mesh(filename, mesh_sm);
  CGAL::draw(mesh_sm, filename);

  // Compute the bounding box
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh_sm);
  std::cout << bbox << std::endl;
  Surface_mesh bbox_pos_sm, bbox_neg_sm;
  split_mesh(bbox, bbox_pos_sm, bbox_neg_sm, 1e-5);
  CGAL::draw(bbox_pos_sm, "Bounding box positive");
  CGAL::draw(bbox_neg_sm, "Bounding box negative");

  auto bbox_sm = bbox_pos_sm;

  // Compute the difference between the bounding box and the mesh
  Surface_mesh complement_sm;
  soup_difference(bbox_sm, mesh_sm, complement_sm);
  auto complement_normals = complement_sm.template add_property_map<Fd, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
  CGAL::Polygon_mesh_processing::compute_face_normals(complement_sm, complement_normals, np);
  merge_coplanar_faces(complement_sm, complement_normals, np);
  retriangulate_faces(complement_sm, complement_normals);
  CGAL::draw(complement_sm, "Complement sm");

  Polyhedral_mesh complement_pm;
  CGAL::copy_face_graph(complement_sm, complement_pm);
  CGAL::draw(complement_pm, "Complement pm");

  // Decompose the difference into convex pieces and compute their Gaussian maps
  std::vector<Arrangement> gms1;
  gaussian_maps(complement_pm, std::back_inserter(gms1), kernel);
  std::cout << "Complement decomposed into " << gms1.size() << " pieces\n";

  // double size = 1e-5;
  double size = 0.2;
  auto cube_iso = make_iso_cube(size, kernel);
  auto cube_pm = to_mesh<Polyhedral_mesh>(cube_iso);

  CGAL::draw(cube_pm, "Iso Cube pm");
  Arrangement cube_gm = gaussian_map<Arrangement>(cube_pm, kernel);
  auto ms_pm = minkowski_sum_3<Polyhedral_mesh>(gms1.begin(), gms1.end(), cube_gm);
  CGAL::draw(ms_pm, "Minkowski Sum pm");

  Surface_mesh ms_sm;
  CGAL::copy_face_graph(ms_pm, ms_sm);

  auto normals = ms_sm.template add_property_map<Fd, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
  CGAL::Polygon_mesh_processing::compute_face_normals(ms_sm, normals, np);
  merge_coplanar_faces(ms_sm, normals, np);
  retriangulate_faces(ms_sm, normals);
  CGAL::draw(ms_sm, "Minkowski Sum sm");

  Surface_mesh result;
  corefine_difference(bbox_sm, ms_sm, result);
  CGAL::draw(result, "result");

  return 0;
}
