#include <vector>

#include <boost/program_options.hpp>

#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/boost/graph/copy_face_graph.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/convex_decomposition_3.h>
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

#include "cgalex/Adder.h"
#include "cgalex/arr_gaussian_map.h"
#include "cgalex/corefine_difference.h"
#include "cgalex/corefine_union.h"
#include "cgalex/Face_normal_map.h"
#include "cgalex/find_file_fullname.h"
#include "cgalex/io_paths.h"
#include "cgalex/gaussian_map.h"
#include "cgalex/merge_coplanar_faces.h"
#include "cgalex/Paths.h"
#include "cgalex/retriangulate_faces.h"
#include "cgalex/soup_difference.h"
#include "cgalex/soup_union.h"

#include "Extended_polyhedron_items.h"
#include "Polyhedron_builder.h"

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

//! \brief obtains the default value of the input path
const Path def_input_path() {
  static const Path s_def_input_path(".");
  return s_def_input_path;
}

//! \brief finds the input file.
std::string find_input_file(const po::variables_map& vm, const std::string& filename)
{ return find_file_fullname(vm["input-path"].as<Paths>(), filename); }

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

template <typename Kernel_>
typename Kernel::Iso_cuboid_3 make_iso_cube(double size, const Kernel_& kernel) {
  using Kernel = Kernel_;

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
  points_pos[1] = Point_3(midx, bbox.ymin()-offset, bbox.zmin()-offset);
  points_pos[2] = Point_3(midx, bbox.ymax()+offset, bbox.zmin()-offset);
  points_pos[3] = Point_3(bbox.xmin()-offset, bbox.ymax()+offset, bbox.zmin()-offset);
  points_pos[4] = Point_3(bbox.xmin()-offset, bbox.ymin()-offset, bbox.zmax()+offset);
  points_pos[5] = Point_3(midx, bbox.ymin()-offset, bbox.zmax()+offset);
  points_pos[6] = Point_3(midx, bbox.ymax()+offset, bbox.zmax()+offset);
  points_pos[7] = Point_3(bbox.xmin()-offset, bbox.ymax()+offset, bbox.zmax()+offset);
  CGAL::convex_hull_3(points_pos.begin(), points_pos.end(), mesh_pos);

  std::vector<Point_3> points_neg(8);
  points_neg[0] = Point_3(midx, bbox.ymin()-offset, bbox.zmin()-offset);
  points_neg[1] = Point_3(bbox.xmax()+offset, bbox.ymin()-offset, bbox.zmin()-offset);
  points_neg[2] = Point_3(bbox.xmax()+offset, bbox.ymax()+offset, bbox.zmin()-offset);
  points_neg[3] = Point_3(midx, bbox.ymax()+offset, bbox.zmin()-offset);
  points_neg[4] = Point_3(midx, bbox.ymin()-offset, bbox.zmax()+offset);
  points_neg[5] = Point_3(bbox.xmax()+offset, bbox.ymin()-offset, bbox.zmax()+offset);
  points_neg[6] = Point_3(bbox.xmax()+offset, bbox.ymax()+offset, bbox.zmax()+offset);
  points_neg[7] = Point_3(midx, bbox.ymax()+offset, bbox.zmax()+offset);
  CGAL::convex_hull_3(points_neg.begin(), points_neg.end(), mesh_neg);

  return true;
}

/*!
 */
template <typename Mesh, typename OutputIterator, typename Kernel_>
OutputIterator gaussian_maps(Mesh& mesh, OutputIterator oi, const Kernel_& kernel) {
  std::cout << "gaussian_maps()\n";
  using Kernel = Kernel_;

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
    // draw(mesh, "Convex part");
    *oi++ = gaussian_map<Arrangement>(mesh, kernel);
  }
  return oi;
}

/*!
 */
template <typename Mesh_, typename InputIterator, typename Arrangement_>
void minkowski_sum_3(InputIterator begin, InputIterator end, const Arrangement_& gm2, Mesh_& ms) {
  std::cout << "minkowski_sum_3()\n";
  using Arrangement = Arrangement_;
  using Mesh = Mesh_;

  std::size_t i = 0;
  CGAL::Arr_face_overlay_traits<Arrangement, Arrangement, Arrangement, Adder> overlay_traits;
  bool first = true;
  for (auto it = begin; it != end; ++it) {
    std::cout << "i: " << i++ << std::endl;
    Arrangement gm;
    CGAL::overlay(*it, gm2, gm, overlay_traits);
    using Halfedge_ds = typename Mesh::HalfedgeDS;
    Polyhedron_builder<Halfedge_ds, Arrangement> surface(gm);
    if (first) {
      ms.delegate(surface);
      PMP::triangulate_faces(ms);
      first = false;
      continue;
    }
    Mesh mesh;
    mesh.delegate(surface);
    PMP::triangulate_faces(mesh);
    corefine_union(ms, mesh, ms);
  }
}

/*!
 */
template <typename Mesh, typename Arrangement_>
void extract_mesh(const Arrangement_& gm, Mesh& mesh) {
  using Arrangement = Arrangement_;
  using Halfedge_ds = typename Mesh::HalfedgeDS;

  Polyhedron_builder<Halfedge_ds, Arrangement> surface(gm);
  mesh.delegate(surface);
  PMP::triangulate_faces(mesh);
}

/*!
 */
template <typename Mesh_, typename Kernel_>
void triangulate(Mesh_& mesh, const Kernel_& kernel) {
  using Mesh = Mesh_;

  auto np = params::geom_traits(kernel);
  Face_normal_map<Mesh> normals;
  CGAL::Polygon_mesh_processing::compute_face_normals(mesh, normals, np);
  merge_coplanar_faces(mesh, normals, np);
  triangulate_faces(mesh, normals);
}

/*!
 */
template <typename Mesh_, typename Arrangement_, typename Kernel_>
void contract_mesh(Mesh_& bbox_part, Mesh_& mesh, const Arrangement_& cube_gm, Mesh_& contracted_part,
                   const Kernel_& kernel) {
  using Mesh = Mesh_;
  using Arrangement = Arrangement_;
  using Kernel = Kernel_;

  auto np = params::geom_traits(kernel);

  // Compute the difference between the bounding box and the mesh
  Mesh complement;
  soup_difference(bbox_part, mesh, complement);

#if USING_SURFACE_MESH
  auto complement_normals = complement_sm.template add_property_map<Fd, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
#else
  Face_normal_map<Mesh> complement_normals;
#endif
  CGAL::Polygon_mesh_processing::compute_face_normals(complement, complement_normals, np);
  merge_coplanar_faces(complement, complement_normals, np);
  triangulate_faces(complement, complement_normals);
  // CGAL::draw(complement, "Complement");

  // Polyhedral_mesh complement;
  // CGAL::copy_face_graph(complement_sm, complement);
  // CGAL::draw(complement, "Complement");

  // Decompose the difference into convex pieces and compute their Gaussian maps
  std::vector<Arrangement> gms1;
  gaussian_maps(complement, std::back_inserter(gms1), kernel);
  std::cout << "Complement decomposed into " << gms1.size() << " pieces\n";

  // Compute the contracted part
#if USING_SURFACE_MESH
  CGAL::copy_face_graph(bbox_part, contracted_part);
#else
  contracted_part = bbox_part;
#endif

  std::size_t i = 0;
  CGAL::Arr_face_overlay_traits<Arrangement, Arrangement, Arrangement, Adder> overlay_traits;
  for (auto it = gms1.begin(); it != gms1.end(); ++it) {
    std::cout << "i: " << i++ << std::endl;
    Arrangement gm;
    CGAL::overlay(*it, cube_gm, gm, overlay_traits);
    Mesh ms;
    extract_mesh(gm, ms);
    corefine_difference(contracted_part, ms, contracted_part);
    PMP::triangulate_faces(contracted_part);
  }
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

  std::size_t verbose = 0;
  auto do_draw = false;
  double size = 0.1;
  std::string input_dir = ".";
  std::string fullname;

  try {
    // Define options
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("verbose,v", po::value<std::size_t>(&verbose)->implicit_value(1), "set verbosity level (0 = quiet")
      ("size,s", po::value<double>(&size)->default_value(0.1), "set cube size")
      ("draw,d", po::value<bool>(&do_draw)->implicit_value(false), "draw the meshes")
      ("input-path,i", po::value<Paths>()->composing()->default_value({def_input_path()}), "input path")
      ("filename", po::value<std::string>(), "First file name")
    ;

    po::positional_options_description p;
    std::string filename;
    p.add("filename", 1);

    // Parse options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
    po::notify(vm);

    // Help
    if (vm.count("help")) {
      std::cout << desc << "\n";
      return 0;
    }

    // Verbosity: show directory being scanned
    if (verbose > 0) {
      std::cout << "Searching directory: " << input_dir << "\n";
    }

    if (! vm.count("filename")) {
      // throw ...
      std::cout << "Missing first filename \n";
      return 1;
    }
    filename = vm["filename"].as<std::string>();
    fullname = find_input_file(vm, filename);
    if (fullname.empty()) {
      throw std::logic_error(std::string("Error: Cannot find file ").append(filename));
      return 1;
    }
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what() << "\n";
    return 1;
  }

  Kernel kernel;

  // Generate an epsilon cube
  // double size = 1e-5;
  auto cube_iso = make_iso_cube(size, kernel);
  auto cube_pm = to_mesh<Polyhedral_mesh>(cube_iso);

  // Compute the Gaussian map of the cube.
  Arrangement cube_gm = gaussian_map<Arrangement>(cube_pm, kernel);

  // Construct the mesh.
  Polyhedral_mesh mesh;
  auto rc1 = CGAL::IO::read_polygon_mesh(fullname, mesh);
  if (do_draw) CGAL::draw(mesh, fullname.c_str());

  // Compute the bounding box and split the mesh.
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  std::cout << bbox << std::endl;
  Polyhedral_mesh bbox_pos, bbox_neg;
  split_mesh(bbox, bbox_pos, bbox_neg, 1e-5);

  // Contract.
  Polyhedral_mesh contracted_pos, contracted_neg, contracted;
  contract_mesh(bbox_pos, mesh, cube_gm, contracted_pos, kernel);
  if (do_draw) CGAL::draw(contracted_pos, "Contracted pos");
  contract_mesh(bbox_neg, mesh, cube_gm, contracted_neg, kernel);
  if (do_draw) CGAL::draw(contracted_pos, "Contracted neg");

  corefine_union(contracted_pos, contracted_neg, contracted);
  triangulate(contracted, kernel);
  if (do_draw) CGAL::draw(contracted, "Contracted");

  // Expand.
  std::vector<Arrangement> gms;
  gaussian_maps(contracted, std::back_inserter(gms), kernel);
  std::cout << "contracted decomposed into " << gms.size() << " pieces\n";
  Polyhedral_mesh expanded;
  minkowski_sum_3<Polyhedral_mesh>(gms.begin(), gms.end(), cube_gm, expanded);
  triangulate(expanded, kernel);
  if (do_draw) CGAL::draw(expanded, "Expanded");
  CGAL::IO::write_polygon_mesh("expanded.off", expanded, CGAL::parameters::stream_precision(17));

  return 0;
}
