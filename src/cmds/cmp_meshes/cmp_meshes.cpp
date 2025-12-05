#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>

#include <boost/program_options.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/graph/isomorphism.hpp>
#include <boost/graph/vf2_sub_graph_iso.hpp>

#include <CGAL/Basic_viewer.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/convexity_check_3.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>

#include "merge_coplanar_faces.h"
#include "triangulate_faces.h"

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
// using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using NT = Kernel::FT;

using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;
using Triangle = boost::container::small_vector<std::size_t, 3>;

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

// Pass indices by value
void apply_permutation(std::vector<Mesh>& v, std::vector<int> indices) {
  for (size_t i = 0; i < indices.size(); i++) {
    auto current = i;
    while (i != indices[current]) {
      auto next = indices[current];
      std::swap(v[current], v[next]);
      indices[current] = current;
      current = next;
    }
    indices[current] = current;
  }
}

//
bool read_mesh(const std::string& filename, Mesh& mesh, const Kernel& kernel) {
  // params::verbose(true).repair_polygon_soup(true)
  if (! PMP::IO::read_polygon_mesh(filename, mesh)) {
    std::cerr << "Error: Cannot read " << filename << "\n";
    return false;
  }

  auto is_valid = mesh.is_valid();
  if (! is_valid) {
    std::cerr << "Error: Mesh " << filename << " is invalid\n";
    return false;
  }

  auto is_closed = CGAL::is_closed(mesh);
  if (! is_closed) {
    std::cerr << "Error: Mesh " << filename << " is not closed\n";
    return false;
  }

  auto is_tri = CGAL::is_triangle_mesh(mesh);
  if (! is_tri) PMP::triangulate_faces(mesh, params::geom_traits(kernel));

  PMP::remove_isolated_vertices(mesh);

  using Vector_3 = typename Kernel::Vector_3;
  auto np = CGAL::parameters::geom_traits(kernel);
  auto normals = mesh.add_property_map<face_descriptor, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
  PMP::compute_face_normals(mesh, normals, np);
  merge_coplanar_faces(mesh, normals, np);

  return true;
}

// Callback used by vf2 to signal success
struct isomorphism_callback {
  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  // true to stop search, false to continue looking for all isomorphisms
  bool operator()(CorrespondenceMap1To2, CorrespondenceMap2To1) { return true; }
};

bool compare(const Mesh& mesh1, const Mesh& mesh2, const Kernel& kernel) {
  // 2. Compare number of features.
  auto num_vertices1 = mesh1.number_of_vertices();
  auto num_edges1 = mesh1.number_of_edges();
  auto num_faces1 = mesh1.number_of_faces();

  auto num_vertices2 = mesh2.number_of_vertices();
  auto num_edges2 = mesh2.number_of_edges();
  auto num_faces2 = mesh2.number_of_faces();

  if (num_vertices1 != num_vertices2) {
    std::cout << "Meshes differ in number of vertices " << num_vertices1 << ", " << num_vertices2 << "\n";
    return false;
  }
  std::cout << "# of vertices: " << num_vertices1 << "\n";

  if (num_edges1 != num_edges2) {
    std::cout << "Meshes differ in number of edges " << num_edges1 << ", " << num_edges2 << "\n";
    return false;
  }
  std::cout << "# of edges: " << num_edges1 << "\n";

  if (num_faces1 != num_faces2) {
    std::cout << "Meshes differ in number of faces " << num_faces1 << ", " << num_faces2 << "\n";
    return false;
  }
  std::cout << "# of faces: " << num_faces1 << "\n";

  // 2. Compare volume and surface area
  auto volume1 = PMP::volume(mesh1);
  auto volume2 = PMP::volume(mesh2);
  if (volume1 != volume2) {
    std::cout << "Meshes differ in volume " << volume1 << ", " << volume2 << "\n";
    return false;
  }
  std::cout << "Volume: " << volume1 << "\n";

  auto area1 = PMP::area(mesh1);
  auto area2 = PMP::area(mesh2);
  if (area1 != area2) {
    std::cout << "Meshes differ in boundary area " << area1 << ", " << area2 << "\n";
    return false;
  }
  std::cout << "Boundary area: " << area1 << "\n";

  // 3. Compare characteristics (i.e., self-intersection, convexity)

  auto self_intersect1 = PMP::does_self_intersect(mesh1);
  auto self_intersect2 = PMP::does_self_intersect(mesh2);
  if (self_intersect1 != self_intersect2) {
    std::cout << "Meshes differ in self intersection " << self_intersect1 << ", " << self_intersect2 << "\n";
    return false;
  }
  std::cout << "Self intersect: " << self_intersect1 << "\n";

  auto is_convex1 = CGAL::is_strongly_convex_3(mesh1, kernel);
  auto is_convex2 = CGAL::is_strongly_convex_3(mesh2, kernel);
  if (is_convex1 != is_convex2) {
    std::cout << "Meshes differ in convexity " << is_convex1 << ", " << is_convex2 << "\n";
    return false;
  }
  std::cout << "Is convex: " << is_convex1 << "\n";

  /* Compare the widths (obb sizes)
   */

  /* Compare isomorphism
   */
  // boost::isomorphism_map<Mesh, Mesh> isom_map;
  bool result = boost::vf2_graph_iso(mesh1, mesh2, isomorphism_callback());
  if (! result) {
    std::cout << "Meshes are not isomorphic\n";
    return false;
  }

  return true;
}

//! Main entry
int main(int argc, char* argv[]) {
  bool do_draw = false;
  std::size_t verbose = 0;
  std::string input_dir = ".";
  std::string filename1;
  std::string filename2;

  try {
    // Define options
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("verbose,v", po::value<std::size_t>(&verbose)->implicit_value(1),
       "set verbosity level (0 = quiet, 1 = show directory, 2 = also count files)")
      ("draw,d", po::value<bool>(&do_draw)->implicit_value(false), "set draw")
      ("input-directory,i", po::value<std::string>(&input_dir)->default_value("."),
       "input directory (default: current directory)")
      ("filename1", po::value<std::string>(), "First file name")
      ("filename2", po::value<std::string>(), "Second file name");
    ;

    po::positional_options_description p;
    p.add("filename1", 1);
    p.add("filename2", 1);

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
    if (verbose >= 1) {
      std::cout << "Scanning directory: " << input_dir << "\n";
    }

    if (! vm.count("filename1")) {
      // throw ...
      std::cerr << "Missing first filename \n";
      return 1;
    }
    filename1 = vm["filename1"].as<std::string>();

    if (! vm.count("filename2")) {
      // throw ...
      std::cerr << "Missing second filename \n";
      return 1;
    }
    filename2 = vm["filename2"].as<std::string>();
  }
  catch (std::exception& e) {
    std::cerr << "Exception: " << e.what() << "\n";
    return 1;
  }

  Kernel kernel;

  Mesh mesh1;
  if (! read_mesh(filename1, mesh1, kernel)) return -1;
  if (mesh1.is_empty()) {
    std::cout << "Mesh " << filename1 << " is empty" << "\n";
    return 0;
  }
  Mesh mesh2;
  if (! read_mesh(filename2, mesh2, kernel)) return -1;
  if (mesh2.is_empty()) {
    std::cout << "Mesh " << filename2 << " is empty" << "\n";
    return 0;
  }

  CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso;
  gso.ignore_all_vertices(true);
  gso.ignore_all_edges(true);
  gso.colored_face = [](const Mesh&, face_descriptor) -> bool { return true; };
  gso.face_color = [] (const Mesh&, face_descriptor fh) -> CGAL::IO::Color {
    if (fh == boost::graph_traits<Mesh>::null_face()) return CGAL::IO::Color(100, 125, 200);
    return get_random_color(CGAL::get_default_random());
  };
  CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso1;
  gso1.ignore_all_vertices(true);
  gso1.ignore_all_edges(true);
  gso1.colored_face = [](const Mesh&, face_descriptor) -> bool { return true; };
  gso1.face_color = [](const Mesh&, face_descriptor fh) -> CGAL::IO::Color { return CGAL::IO::Color(0, 0, 255); };
  CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso2;
  gso2.ignore_all_vertices(true);
  gso2.ignore_all_edges(true);
  gso2.colored_face = [](const Mesh&, face_descriptor) -> bool { return true; };
  gso2.face_color = [](const Mesh&, face_descriptor fh) -> CGAL::IO::Color { return CGAL::IO::Color(255, 0, 0); };
  CGAL::Graphics_scene scene;
  CGAL::add_to_graphics_scene(mesh1, scene, gso1);
  CGAL::add_to_graphics_scene(mesh2, scene, gso2);
  CGAL::draw_graphics_scene(scene);

  // 1. Compare # of components

  auto fccmap1 = mesh1.add_property_map<face_descriptor, std::size_t>("f:CC").first;
  auto num_ccs1 = PMP::connected_components(mesh1, fccmap1);
  auto fccmap2 = mesh2.add_property_map<face_descriptor, std::size_t>("f:CC").first;
  auto num_ccs2 = PMP::connected_components(mesh2, fccmap2);
  if (num_ccs1 != num_ccs2) {
    std::cout << "Meshes differ in # of components " << num_ccs1 << ", " << num_ccs2 << "\n";
    return 0;
  }
  std::cout << "# of components: " << num_ccs1 << "\n";

  // Split connected components
  std::vector<Mesh> meshes1;
  meshes1.reserve(num_ccs1);
  PMP::split_connected_components(mesh1, meshes1);
  std::vector<Mesh> meshes2;
  meshes2.reserve(num_ccs2);
  PMP::split_connected_components(mesh2, meshes2);

  // Compute polygon soup and sort
  using Points = std::vector<Kernel::Point_3>;

  std::vector<Points> points_of_meshes1(num_ccs1);
  for (std::size_t i = 0; i < num_ccs1; ++i) {
    auto& mesh = meshes1[i];
    PMP::remove_isolated_vertices(mesh);
    points_of_meshes1[i].reserve(mesh.number_of_vertices());
    for (const auto& p : mesh.points()) {
      points_of_meshes1[i].push_back(p);
      // std::cout << p << " ";
    }
    // std::cout << std::endl;
    std::sort(points_of_meshes1[i].begin(), points_of_meshes1[i].end());
  }
  std::vector<int> indices1(num_ccs1);
  std::iota(indices1.begin(), indices1.end(), 0);
  std::sort(indices1.begin(), indices1.end(), [&](int i, int j){ return points_of_meshes1[i] < points_of_meshes1[j]; });
  apply_permutation(meshes1, indices1);

  std::vector<Points> points_of_meshes2(num_ccs2);
  for (std::size_t i = 0; i < num_ccs2; ++i) {
    auto& mesh = meshes2[i];
    // PMP::remove_isolated_vertices(mesh);
    points_of_meshes2[i].reserve(mesh.number_of_vertices());
    for (const auto& p : mesh.points()) {
      points_of_meshes2[i].push_back(p);
      // std::cout << p << " ";
    }
    // std::cout << std::endl;
    std::sort(points_of_meshes2[i].begin(), points_of_meshes2[i].end());
  }
  std::vector<int> indices2(num_ccs2);
  std::iota(indices2.begin(), indices2.end(), 0);
  std::sort(indices2.begin(), indices2.end(), [&](int i, int j){ return points_of_meshes2[i] < points_of_meshes2[j]; });
  apply_permutation(meshes2, indices2);

  // Compare meshes
  for (std::size_t i = 0; i < num_ccs1; ++i) {
    const auto& points1 = points_of_meshes1[indices1[i]];
    const auto& points2 = points_of_meshes2[indices2[i]];
    if (points1.size() != points2.size()) {
      std::cout << "Sub meshes " << indices1[i] << " and " << indices2[i] << " differ in number of points "
                << points1.size() << ", " << points2.size() << "\n";
      return false;
    }
    if (points1 != points2) {
      for (const auto& p : points1) std::cout << p << " ";
      std::cout << std::endl;
      for (const auto& p : points2) std::cout << p << " ";
      std::cout << std::endl;
      std::cout << "Sub meshes " << indices1[i] << " and " << indices2[i] << " differ in points\n";
      return false;
    }
    const auto& mesh1 = meshes1[i];
    const auto& mesh2 = meshes2[i];
    auto diff = compare(mesh1, mesh2, kernel);
  }

  /* Use Side_of_triangle_mesh.
   * To begin with, you might want to compute the volume and the area.
   * You can also use PMP::approximate_Hausdorff_distance and
   * bounded_error_Hausdorff_distance()
   */

  std::cout << "The meshes are isomorphic\n";
  return 0;
}
