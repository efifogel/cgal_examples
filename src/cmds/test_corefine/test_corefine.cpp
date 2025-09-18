#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <array>

#include <boost/program_options.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/IO/polygon_soup_io.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
// using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;

using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

struct Vector_pmap_wrapper {
  std::vector<bool>& vect;
  Vector_pmap_wrapper(std::vector<bool>& v) : vect(v) {}
  friend bool get(const Vector_pmap_wrapper& m, face_descriptor f) { return m.vect[f]; }
  friend void put(const Vector_pmap_wrapper& m, face_descriptor f, bool b) { m.vect[f] = b; }
};

// enum class Op { ALL = 0, DIFFERENCE, INTERSECTION, UNION };
using Op = PMP::Corefinement::Boolean_operation_type;

//! Names of edge operations.
static const char* s_op_names[] = { "union", "intersection", "tm1_minus_tm2", "tm2_minus_tm1", "none" };

//! Convert Op enumeration to std::size_t
inline std::size_t to_index(Op op) { return static_cast<std::size_t>(op); }

namespace std {

template <typename OutputStream>
inline OutputStream& operator<<(OutputStream& os, Op op) {
  os << s_op_names[to_index(op)];
  return os;
}

//! Import a thread type from an input stream.
template <typename InputStream>
inline InputStream& operator>>(InputStream& is, Op& op) {
  std::string name;
  is >> name;
  auto it = std::find(std::begin(s_op_names), std::end(s_op_names), name);
  if (it == std::end(s_op_names)) throw std::invalid_argument(name);
  auto diff = std::distance(std::begin(s_op_names), it);
  op = static_cast<Op>(diff);
  return is;
}

}

//! Main entry
int main(int argc, char* argv[]) {
  bool do_draw = false;
  std::size_t verbose = 0;
  std::string input_dir = ".";
  std::string filename1;
  std::string filename2;
  bool corefine = false;
  Op op = PMP::Corefinement::TM1_MINUS_TM2;

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
      ("op,o", po::value<Op>(&op)->default_value(PMP::Corefinement::TM1_MINUS_TM2), "operation")
      ("corefine,c", po::value<bool>(&corefine)->implicit_value(false), "set corefine")
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

  Mesh mesh1, mesh2;
  if (! PMP::IO::read_polygon_mesh(filename1, mesh1) || ! PMP::IO::read_polygon_mesh(filename2, mesh2)) {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso;
  gso.ignore_all_vertices(true);
  gso.ignore_all_edges(true);
  gso.colored_face = [](const Mesh&, typename boost::graph_traits<Mesh>::face_descriptor) -> bool
  { return true; };
  gso.face_color =  [] (const Mesh&, typename boost::graph_traits<Mesh>::face_descriptor fh) -> CGAL::IO::Color {
    if (fh == boost::graph_traits<Mesh>::null_face()) return CGAL::IO::Color(100, 125, 200);
    return get_random_color(CGAL::get_default_random());
  };
  CGAL::Graphics_scene scene;
  CGAL::add_to_graphics_scene(mesh1, scene, gso);
  CGAL::add_to_graphics_scene(mesh2, scene, gso);
  CGAL::draw_graphics_scene(scene);

  // create a property on edges to indicate whether they are constrained
  Mesh::Property_map<edge_descriptor,bool> is_constrained_map =
    mesh1.add_property_map<edge_descriptor,bool>("e:is_constrained", false).first;

  PMP::Corefinement::Non_manifold_output_visitor<Mesh> visitor(mesh1, mesh2);
  Mesh result;
  bool valid_result;

  if (corefine) {
    std::array<std::optional<Mesh*>, 4> results;
    std::array<bool, 4> valid_results;
    switch (op) {
     case PMP::Corefinement::UNION:
      results[PMP::Corefinement::INTERSECTION] = &result;
      valid_results = PMP::corefine_and_compute_boolean_operations(mesh1, mesh2, results,
                                                                   params::default_values(), // mesh1 np
                                                                   params::default_values(), // mesh2 np
                                                                   std::make_tuple(params::visitor(visitor).vertex_point_map(get(boost::vertex_point, result)),
                                                                                   params::default_values(),
                                                                                   params::default_values(),
                                                                                   params::default_values()));
      break;

     case PMP::Corefinement::INTERSECTION:
      results[PMP::Corefinement::UNION] = &result;
      valid_results = PMP::corefine_and_compute_boolean_operations(mesh1, mesh2, results,
                                                                   params::default_values(), // mesh1 np
                                                                   params::default_values(), // mesh2 np
                                                                   std::make_tuple(params::default_values(),
                                                                                   params::visitor(visitor).vertex_point_map(get(boost::vertex_point, result)),
                                                                                   params::default_values(),
                                                                                   params::default_values()));
      break;

     case PMP::Corefinement::TM1_MINUS_TM2:
      results[PMP::Corefinement::TM1_MINUS_TM2] = &result;
      valid_results = PMP::corefine_and_compute_boolean_operations(mesh1, mesh2, results,
                                                                   params::default_values(), // mesh1 np
                                                                   params::default_values(), // mesh2 np
                                                                   std::make_tuple(params::default_values(),
                                                                                   params::default_values(),
                                                                                   params::visitor(visitor).vertex_point_map(get(boost::vertex_point, result)),
                                                                                   params::default_values()));
      break;

     case PMP::Corefinement::TM2_MINUS_TM1:
      results[PMP::Corefinement::TM2_MINUS_TM1] = &result;
      valid_results = PMP::corefine_and_compute_boolean_operations(mesh1, mesh2, results,
                                                                   params::default_values(), // mesh1 np
                                                                   params::default_values(), // mesh2 np
                                                                   std::make_tuple(params::default_values(),
                                                                                   params::default_values(),
                                                                                   params::default_values(),
                                                                                   params::visitor(visitor).vertex_point_map(get(boost::vertex_point, result))));
      break;
    }

    valid_result = valid_results[op];
  }
  else {
    switch (op) {
     case PMP::Corefinement::UNION:
      valid_result = PMP::corefine_and_compute_union(mesh1, mesh2, result, params::visitor(visitor));
      break;

     case PMP::Corefinement::INTERSECTION:
      valid_result = PMP::corefine_and_compute_intersection(mesh1, mesh2, result, params::visitor(visitor));
      break;

     case PMP::Corefinement::TM1_MINUS_TM2:
      valid_result = PMP::corefine_and_compute_difference(mesh1, mesh2, result, params::visitor(visitor));
      break;

     case PMP::Corefinement::TM2_MINUS_TM1:
      valid_result = PMP::corefine_and_compute_difference(mesh2, mesh1, result, params::visitor(visitor));
      break;
    }
  }

  if (valid_result) {
    std::cout << "Difference was successfully computed\n";
    auto self_intersect = PMP::duplicate_non_manifold_vertices(result);
    if (self_intersect) {
      std::cout << "The output self intersect as a result of duplicating non-manifold vertices\n";
    }
  }
  else {
    std::cout << "Difference failed; resorting to handling a non-manifild result\n";
    std::vector<Kernel::Point_3> points;
    std::vector<std::array<std::size_t, 3>> polygons;

    visitor.extract_tm1_minus_tm2(points, polygons);
    // CGAL::IO::write_polygon_soup("inter_soup.off", points, polygons, CGAL::parameters::stream_precision(17));
    // make the soup topologically manifold (but geometrically self-intersecting)
    PMP::orient_polygon_soup(points, polygons);
    // fill a mesh with the intersection
    PMP::polygon_soup_to_polygon_mesh(points, polygons, result);
  }

  CGAL::draw(result, "result");
  CGAL::IO::write_polygon_mesh("result.off", result, CGAL::parameters::stream_precision(17));

  return 0;
}
