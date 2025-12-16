#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <string>
#include <vector>
#include <cmath>

#include <boost/program_options.hpp>
#include <boost/container/small_vector.hpp>
// #include <boost/graph/vf2_sub_graph_iso.hpp>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/convexity_check_3.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/repair_degeneracies.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>

#if defined(CGALEX_HAS_VISUAL)
#include <CGAL/Basic_viewer.h>
#include <CGAL/draw_surface_mesh.h>
#endif

#include "cgalex/find_file_fullname.h"
#include "cgalex/io_paths.h"
#include "cgalex/merge_coplanar_faces.h"
#include "cgalex/Paths.h"
#include "cgalex/triangulate_faces.h"

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
bool read_mesh(const std::string& filename, Mesh& mesh, bool do_repair, std::size_t verbose, const Kernel& kernel) {
  // params::verbose(true).repair_polygon_soup(true)
  if (! PMP::IO::read_polygon_mesh(filename, mesh, params::verbose(true).repair_polygon_soup(true))) {
    std::cout << "Error: Cannot read " << filename << "\n";
    return false;
  }

  auto is_valid = mesh.is_valid();
  if (! is_valid) {
    std::cout << "Error: Mesh " << filename << " is invalid\n";
    return false;
  }

  auto is_closed = CGAL::is_closed(mesh);
  if (! is_closed) {
    std::cout << "Error: Mesh " << filename << " is not closed\n";
    return false;
  }

  auto is_tri = CGAL::is_triangle_mesh(mesh);
  if (! is_tri) PMP::triangulate_faces(mesh, params::geom_traits(kernel));

  if (do_repair) {
    if (verbose > 0) std::cout << "Repairing\n";
    PMP::remove_degenerate_faces(mesh);
    PMP::remove_degenerate_edges(mesh);
    PMP::remove_isolated_vertices(mesh);
  }

  return true;
}

// Callback used by vf2 to signal success
struct isomorphism_callback {
  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  // true to stop search, false to continue looking for all isomorphisms
  bool operator()(CorrespondenceMap1To2, CorrespondenceMap2To1) { return true; }
};

enum class Criterion : int {
  NUMBER_OF_POINTS = 1,
  POINTS,
  NUMBER_OF_VERTICES,
  NUMBER_OF_EDGES,
  NUMBER_OF_FACES,
  VOLUME,
  AREA,
  SELF_INTERSECTION,
  CONVEXITY,
  ISOMORPHISM
};

struct Result {
  bool m_ok = true;
  Criterion m_criterion;
  std::string m_message;
  Result() : m_ok(true) {}
  Result(Criterion criterion, const std::string& message) : m_ok(false), m_criterion(criterion), m_message(message) {}
  Result(bool ok, Criterion criterion, const std::string& message) : m_ok(ok), m_criterion(criterion), m_message(message) {}
};

//
std::vector<Result> compare(const Mesh& mesh1, const Mesh& mesh2, const Kernel& kernel, double tolerance, std::size_t verbose) {
  std::vector<Result> results;

  // 1. Compare number of features.
  auto num_vertices1 = mesh1.number_of_vertices();
  auto num_vertices2 = mesh2.number_of_vertices();
  if (num_vertices1 != num_vertices2) {
    results.emplace_back(Criterion::NUMBER_OF_VERTICES, std::to_string(num_vertices1) + ", " + std::to_string(num_vertices2));
  }
  else if (verbose > 1) results.emplace_back(true, Criterion::NUMBER_OF_VERTICES, std::to_string(num_vertices1));

  auto num_edges1 = mesh1.number_of_edges();
  auto num_edges2 = mesh2.number_of_edges();
  if (num_edges1 != num_edges2) {
    results.emplace_back(Criterion::NUMBER_OF_EDGES, std::to_string(num_edges1) + ", " + std::to_string(num_edges2));
  }
  else if (verbose > 1) results.emplace_back(true, Criterion::NUMBER_OF_EDGES, std::to_string(num_edges1));

  auto num_faces1 = mesh1.number_of_faces();
  auto num_faces2 = mesh2.number_of_faces();

  if (num_faces1 != num_faces2) {
    results.emplace_back(Criterion::NUMBER_OF_FACES, std::to_string(num_faces1) + ", " + std::to_string(num_faces2));
  }
  else if (verbose > 1) results.emplace_back(true, Criterion::NUMBER_OF_FACES, std::to_string(num_faces1));

  auto tmesh1 = mesh1;
  using Vector_3 = typename Kernel::Vector_3;
  std::optional<Mesh::Property_map<face_descriptor, Vector_3>> optional_map1 = mesh1.property_map<face_descriptor, Vector_3>("f:normals");
  CGAL_assertion(optional_map1);
  Mesh::Property_map<face_descriptor, Vector_3> normals1 = *optional_map1;
  triangulate_faces(tmesh1, normals1, params::geom_traits(kernel));

  auto tmesh2 = mesh2;
  using Vector_3 = typename Kernel::Vector_3;
  std::optional<Mesh::Property_map<face_descriptor, Vector_3>> optional_map2 = mesh2.property_map<face_descriptor, Vector_3>("f:normals");
  CGAL_assertion(optional_map2);
  Mesh::Property_map<face_descriptor, Vector_3> normals2 = *optional_map2;
  triangulate_faces(tmesh2, normals2, params::geom_traits(kernel));

  // 2. Compare volume and surface area
  auto volume1 = CGAL::to_double(PMP::volume(tmesh1));
  auto volume2 = CGAL::to_double(PMP::volume(tmesh2));
  if (std::abs(volume1 - volume2) > tolerance) {
    results.emplace_back(Criterion::VOLUME, std::to_string(volume1) + ", " + std::to_string(volume2));
  }
  else if (verbose > 1) results.emplace_back(true, Criterion::VOLUME, std::to_string(volume1));

  // Replace with squared_area() of face
  auto area1 = CGAL::to_double(PMP::area(tmesh1));
  auto area2 = CGAL::to_double(PMP::area(tmesh2));
  if (std::abs(area1 - area2) > tolerance) {
    results.emplace_back(Criterion::AREA, std::to_string(area1) + ", " + std::to_string(area2));
  }
  else if (verbose > 1) results.emplace_back(true, Criterion::AREA, std::to_string(area1));

  // 3. Compare characteristics (i.e., self-intersection, convexity)

  auto self_intersect1 = PMP::does_self_intersect(mesh1);
  auto self_intersect2 = PMP::does_self_intersect(mesh2);
  if (self_intersect1 != self_intersect2) {
    results.emplace_back(Criterion::SELF_INTERSECTION, std::to_string(self_intersect1) + ", " + std::to_string(self_intersect2));
  }
  else if (verbose > 1) results.emplace_back(true, Criterion::SELF_INTERSECTION, std::to_string(self_intersect1));

  auto is_convex1 = CGAL::is_strongly_convex_3(mesh1, kernel);
  auto is_convex2 = CGAL::is_strongly_convex_3(mesh2, kernel);
  if (is_convex1 != is_convex2) {
    results.emplace_back(Criterion::CONVEXITY, std::to_string(is_convex1) + ", " + std::to_string(is_convex2));
  }
  else if (verbose > 1) results.emplace_back(true, Criterion::CONVEXITY, std::to_string(is_convex1));

  /* Compare the widths (obb sizes)
   */

  /* Compare isomorphism
   */
  // boost::isomorphism_map<Mesh, Mesh> isom_map;
# if 0
  bool isomorphic = boost::vf2_graph_iso(mesh1, mesh2, isomorphism_callback());
  if (! isomorphic) {
    results.emplace_back(Criterion::ISOMORPHISM, "Unknown Information");
  }
#endif

  std::vector<std::pair<face_descriptor, face_descriptor>> common;
  std::vector<face_descriptor> m1_only, m2_only;
  PMP::match_faces(mesh1, mesh2, std::back_inserter(common), std::back_inserter(m1_only), std::back_inserter(m2_only));
  if (! m1_only.empty() || ! m2_only.empty()) {
    results.emplace_back(Criterion::ISOMORPHISM,
                         std::to_string(m1_only.size()) + " only in 1st, " + std::to_string(m2_only.size()) + " only in 2nd");

    if (verbose > 2) {
      auto vpm1 = get_const_property_map(boost::vertex_point, tmesh1);
      auto vpm2 = get_const_property_map(boost::vertex_point, tmesh2);
      for (auto fd : m1_only) {
        std::cout << fd << "\n";
        auto rep_hd = halfedge(fd, tmesh1);
        for (auto hd : CGAL::halfedges_around_face(rep_hd, tmesh1)) {
          auto vd = CGAL::target(hd, tmesh1);
          const auto& p = get(vpm1, vd);
          std::cout << p << std::endl;
        }
        std::cout << std::endl;
      }
      for (auto fd : m2_only) {
        std::cout << fd << "\n";
        auto rep_hd = halfedge(fd, tmesh2);
        for (auto hd : CGAL::halfedges_around_face(rep_hd, tmesh2)) {
          auto vd = CGAL::target(hd, tmesh2);
          const auto& p = get(vpm2, vd);
          std::cout << p << std::endl;
        }
        std::cout << std::endl;
      }
    }
  }

  return results;
}

//! \brief obtains the default value of the input path
const Path def_input_path() {
  static const Path s_def_input_path(".");
  return s_def_input_path;
}

//! \brief finds the input file.
std::string find_input_file(const po::variables_map& vm, const std::string& filename)
{ return find_file_fullname(vm["input-path"].as<Paths>(), filename); }

//! Main entry
int main(int argc, char* argv[]) {
  bool do_draw = false;
  bool do_repair = true;
  std::size_t verbose = 0;
  std::string input_dir = ".";
  double tolerance = 1e-5;
  bool exhaustive = false;
  std::string fullname1;
  std::string fullname2;

  try {
    // Define options
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("verbose,v", po::value<std::size_t>(&verbose)->implicit_value(1), "set verbosity level (0 = quiet")
      ("draw,d", po::value<bool>(&do_draw)->implicit_value(true), "draw the meshes")
      ("exhaustive,e", po::value<bool>(&exhaustive)->implicit_value(true), "exhaustive")
      ("repair,r", po::value<bool>(&do_repair)->implicit_value(true), "repair the meshes")
      ("input-path,i", po::value<Paths>()->composing()->default_value({def_input_path()}), "input path")
      ("tolerance", po::value<double>(&tolerance), "the tolerance when comparing inexact values")
      ("filename1", po::value<std::string>(), "First file name")
      ("filename2", po::value<std::string>(), "Second file name");
    ;

    po::positional_options_description p;
    std::string filename1;
    std::string filename2;
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
    if (verbose > 0) {
      std::cout << "Comparing directory: " << input_dir << "\n";
    }

    if (! vm.count("filename1")) {
      // throw ...
      std::cout << "Missing first filename \n";
      return 1;
    }
    filename1 = vm["filename1"].as<std::string>();
    fullname1 = find_input_file(vm, filename1);
    if (fullname1.empty()) {
      throw std::logic_error(std::string("Error: Cannot find file ").append(filename1));
      return 1;
    }

    if (! vm.count("filename2")) {
      // throw ...
      std::cout << "Error: Missing second filename \n";
      return 1;
    }
    filename2 = vm["filename2"].as<std::string>();
    fullname2 = find_input_file(vm, filename2);
    if (fullname2.empty()) {
      throw std::logic_error(std::string("Error: Cannot find file ").append(filename2));
      return 1;
    }
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what() << "\n";
    return 1;
  }

  Kernel kernel;

  Mesh mesh1, mesh2;
  if (! read_mesh(fullname1, mesh1, do_repair, verbose, kernel)) return -1;
  if (! read_mesh(fullname2, mesh2, do_repair, verbose, kernel)) return -1;

  if (mesh1.is_empty()) {
    if (mesh2.is_empty()) {
      std::cout << "Both meshes are empty\n";
      return 0;
    }
    std::cout << "Mesh " << fullname1 << " is empty" << "\n";
    return -1;
  }
  if (mesh2.is_empty()) {
    std::cout << "Mesh " << fullname2 << " is empty" << "\n";
    return -1;
  }

#if defined(CGALEX_HAS_VISUAL)
  if (do_draw) {
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
  }
#endif

  // 1. Compare # of components

  auto fccmap1 = mesh1.add_property_map<face_descriptor, std::size_t>("f:CC").first;
  auto num_ccs1 = PMP::connected_components(mesh1, fccmap1);
  auto fccmap2 = mesh2.add_property_map<face_descriptor, std::size_t>("f:CC").first;
  auto num_ccs2 = PMP::connected_components(mesh2, fccmap2);
  if (num_ccs1 != num_ccs2) {
    std::cout << "Meshes differ in # of components " << num_ccs1 << ", " << num_ccs2 << "\n";
    return -1;
  }
  if (verbose > 0) std::cout << "# of components: " << num_ccs1 << "\n";

  // Split connected components
  std::vector<Mesh> meshes1;
  meshes1.reserve(num_ccs1);
  PMP::split_connected_components(mesh1, meshes1);
  std::vector<Mesh> meshes2;
  meshes2.reserve(num_ccs2);
  PMP::split_connected_components(mesh2, meshes2);

  // Merge coplanar faces
  using Vector_3 = typename Kernel::Vector_3;
  auto np = CGAL::parameters::geom_traits(kernel);
  for (std::size_t i = 0; i < num_ccs1; ++i) {
    auto& mesh = meshes1[i];
    auto normals = mesh.add_property_map<face_descriptor, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_face_normals(mesh, normals, np);
    merge_coplanar_faces(mesh, normals, np);
  }
  for (std::size_t i = 0; i < num_ccs2; ++i) {
    auto& mesh = meshes2[i];
    auto normals = mesh.add_property_map<face_descriptor, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_face_normals(mesh, normals, np);
    merge_coplanar_faces(mesh, normals, np);
  }

  // Compute polygon soup and sort
  using Points = std::vector<Kernel::Point_3>;

  std::vector<Points> points_of_meshes1(num_ccs1);
  for (std::size_t i = 0; i < num_ccs1; ++i) {
    auto& mesh = meshes1[i];
    PMP::remove_isolated_vertices(mesh);
    points_of_meshes1[i].reserve(mesh.number_of_vertices());
    for (const auto& p : mesh.points()) {
      points_of_meshes1[i].push_back(p);
    }
    // std::cout << std::endl;
    std::sort(points_of_meshes1[i].begin(), points_of_meshes1[i].end());
  }
  std::vector<int> indices1(num_ccs1);
  std::iota(indices1.begin(), indices1.end(), 0);
  std::sort(indices1.begin(), indices1.end(), [&](int i, int j){ return points_of_meshes1[i] < points_of_meshes1[j]; });
  apply_permutation(meshes1, indices1);

  std::vector<Points> points_of_meshes2(num_ccs2);
  if (verbose == 0) std::cout << "Comparing ";
  for (std::size_t i = 0; i < num_ccs2; ++i) {
    if (verbose == 0) std::cout << ".";
    auto& mesh = meshes2[i];
    PMP::remove_isolated_vertices(mesh);
    points_of_meshes2[i].reserve(mesh.number_of_vertices());
    for (const auto& p : mesh.points()) {
      points_of_meshes2[i].push_back(p);
    }
    // std::cout << std::endl;
    std::sort(points_of_meshes2[i].begin(), points_of_meshes2[i].end());
  }
  std::vector<int> indices2(num_ccs2);
  std::iota(indices2.begin(), indices2.end(), 0);
  std::sort(indices2.begin(), indices2.end(), [&](int i, int j){ return points_of_meshes2[i] < points_of_meshes2[j]; });
  apply_permutation(meshes2, indices2);

  // Compare meshes
  bool identical = true;
  std::vector<std::vector<Result>> results(num_ccs1);

  for (std::size_t i = 0; i < num_ccs1; ++i) {
    const auto& points1 = points_of_meshes1[indices1[i]];
    const auto& points2 = points_of_meshes2[indices2[i]];
    if (points1.size() != points2.size()) {
      results[i].emplace_back(Criterion::NUMBER_OF_POINTS, std::to_string(points1.size()) + ", " + std::to_string(points2.size()));
      identical = false;
      if (! exhaustive) break;
    }
    else if (verbose > 1) results[i].emplace_back(true, Criterion::NUMBER_OF_POINTS, std::to_string(points1.size()));

    if (points1 != points2) {
      // if (verbose > 2) {
      //   for (const auto& p : points1) std::cout << p << " ";
      //   std::cout << std::endl;
      //   for (const auto& p : points2) std::cout << p << " ";
      //   std::cout << std::endl;
      // }
      results[i].emplace_back(Criterion::POINTS, "");
      identical = false;
      if (! exhaustive) break;
    }
    else if (verbose > 1) results[i].emplace_back(true, Criterion::POINTS, "");

    const auto& mesh1 = meshes1[i];
    const auto& mesh2 = meshes2[i];

    auto cc_results = compare(mesh1, mesh2, kernel, tolerance, verbose);
    results[i].insert(results[i].end(), cc_results.begin(), cc_results.end());
    if (! results[i].empty()) {
      for (const auto& res : results[i]) {
        if (! res.m_ok) {
          identical = false;
          break;
        }
      }
    }
  }
  if (verbose == 0) std::cout << " " << ((identical) ? "identical" : "different") << std::endl;

  for (std::size_t i = 0; i < num_ccs1; ++i) {
    if (results[i].empty()) continue;
    for (const auto& res : results[i]) {
      if (num_ccs1 == 1) std::cout << "Meshes";
      else std::cout << "Sub meshes " << indices1[i] << " and " << indices2[i];
      if (! res.m_ok) std::cout << " differ in";
      switch (res.m_criterion) {
       case Criterion::NUMBER_OF_POINTS:
        std::cout << " number of points: " << res.m_message << std::endl;
        break;

       case Criterion::POINTS:
        std::cout << " points: " << res.m_message << std::endl;
        break;

       case Criterion::NUMBER_OF_VERTICES:
        std::cout << " number of vertices: " << res.m_message << std::endl;
        break;

       case Criterion::NUMBER_OF_EDGES:
        std::cout << " number of edges: " << res.m_message << std::endl;
        break;

       case Criterion::NUMBER_OF_FACES:
        std::cout << " number of faces: " << res.m_message << std::endl;
        break;

       case Criterion::VOLUME:
        std::cout << " volume: " << res.m_message << std::endl;
        break;

       case Criterion::AREA:
        std::cout << " area: " << res.m_message << std::endl;
        break;

       case Criterion::SELF_INTERSECTION:
        std::cout << " self intersection: " << res.m_message << std::endl;
        break;

       case Criterion::CONVEXITY:
        std::cout << " convexity: " << res.m_message << std::endl;
        break;

       case Criterion::ISOMORPHISM:
        std::cout << " isomorphism: " << res.m_message << std::endl;
        break;

       default: break;
      }
    }
  }
  if (! identical) return -1;

  /* Use Side_of_triangle_mesh.
   * Use PMP::approximate_Hausdorff_distance and bounded_error_Hausdorff_distance()
   */

  // std::cout << "The meshes are isomorphic\n";
  return 0;
}
