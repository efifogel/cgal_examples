#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <list>

#include <boost/program_options.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/draw_arrangement_2.h>
#include <CGAL/Graphics_scene_options.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arr_segment_traits_2<Kernel>;
using Name = std::string;
using Names = std::list<std::string>;
using Dcel = CGAL::Arr_face_extended_dcel<Traits, Names>;
using Arrangement = CGAL::Arrangement_2<Traits, Dcel>;
using Point_2 = Traits::Point_2;
using X_monotone_curve_2 = Traits::X_monotone_curve_2;

using Vertex_const_handle = Arrangement::Vertex_const_handle;
using Halfedge_const_handle = Arrangement::Halfedge_const_handle;
using Face_const_handle = Arrangement::Face_const_handle;

using Vertex_handle = Arrangement::Vertex_handle;
using Halfedge_handle = Arrangement::Halfedge_handle;
using Face_handle = Arrangement::Face_handle;

namespace po = boost::program_options;
namespace fs = std::filesystem;

struct Overlay_traits {
  using V_const_handle = typename Arrangement::Vertex_const_handle;
  using H_const_handle = typename Arrangement::Halfedge_const_handle;
  using F_const_handle = typename Arrangement::Face_const_handle;
  using V_handle = typename Arrangement::Vertex_handle;
  using H_handle = typename Arrangement::Halfedge_handle;
  using F_handle = typename Arrangement::Face_handle;

  void create_face(F_const_handle f1, F_const_handle f2, F_handle f) const {
    f->set_data(f1->data());
    f->data().insert(f->data().end(), f2->data().begin(), f2->data().end());
  }

  void create_vertex(H_const_handle h1, H_const_handle h2, V_handle v) const {}
  void create_vertex(V_const_handle v1, V_const_handle v2, V_handle v) const {}
  void create_vertex(V_const_handle v1, H_const_handle h2, V_handle v) const {}
  void create_vertex(H_const_handle h1, V_const_handle v2, V_handle v) const {}
  void create_vertex(F_const_handle f1, V_const_handle v2, V_handle v) const {}
  void create_vertex(V_const_handle v1, F_const_handle f2, V_handle v) const {}
  void create_edge(H_const_handle h1, H_const_handle h2, H_handle h) const {}
  void create_edge(H_const_handle h1, F_const_handle f2, H_handle h) const {}
  void create_edge(F_const_handle f1, H_const_handle h2, H_handle h) const {}
};

//! Insert a polygon into the arrangement
template <typename InputIterator>
Arrangement polygon_arrangement(const std::string& name, InputIterator begin, InputIterator end) {
  Arrangement arr;
  std::cout << "Adding polygon " << name << "\n";
  auto next = begin;
  auto it = next++;
  for (; next != end; it = next++) {
    X_monotone_curve_2 seg(*it, *next);
    std::cout << "  " << seg << "\n";
    CGAL::insert(arr, seg);
  }
  X_monotone_curve_2 seg(*it, *begin);
  std::cout << "  " << seg << "\n";
  CGAL::insert(arr, seg);
  CGAL_assertion(arr.number_of_faces() == 2);
  auto unbounded_face = arr.unbounded_face();
  CGAL_assertion(unbounded_face->number_of_inner_ccbs() == 1);
  auto inner_ccb_it = unbounded_face->inner_ccbs_begin();
  auto curr = *inner_ccb_it;
  auto inner_face = curr->twin()->face();
  inner_face->data().push_back(name);
  return arr;
}

//! Load all polygons from an input file
bool load_polygons_from_file(const std::string& filename, Arrangement& arr) {
  std::ifstream infile(filename);
  if (! infile.is_open()) {
    std::cerr << "Error: cannot open file " << filename << std::endl;
    return false;
  }

  std::size_t num_polygons;
  infile >> num_polygons;
  if (! infile) {
    std::cerr << "Error: failed to read number of polygons" << std::endl;
    return false;
  }

  for (std::size_t i = 0; i < num_polygons; ++i) {
    std::string name;
    std::size_t num_points;
    infile >> name >> num_points;
    if (! infile) {
      std::cerr << "Error: malformed input while reading polygon header" << std::endl;
      return false;
    }

    std::vector<Point_2> points;
    points.reserve(num_points);

    for (std::size_t j = 0; j < num_points; ++j) {
      double x, y;
      infile >> x >> y;
      if (! infile) {
        std::cerr << "Error: malformed input while reading point "
                  << j << " of polygon " << name << std::endl;
        return false;
      }
      points.emplace_back(x, y);
    }

    Arrangement pgn_arr = polygon_arrangement(name, points.begin(), points.end());
    Arrangement new_arr;
    Overlay_traits overlay_traits;
    CGAL::overlay(arr, pgn_arr, new_arr, overlay_traits);
    arr = new_arr;
  }

  return true;
}

//! The main entry
int main(int argc, char* argv[]) {
  bool do_draw = false;
  std::size_t verbose = 0;
  std::string input_dir = ".";
  std::string filename;

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
      ("filename", po::value<std::string>(), "input file name")
    ;

    po::positional_options_description p;
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
    if (verbose >= 1) {
      std::cout << "Scanning directory: " << input_dir << "\n";
    }

    if (! vm.count("filename")) {
      // throw ...
      std::cerr << "Missing filename \n";
      return 1;
    }
    filename = vm["filename"].as<std::string>();
  }
  catch (std::exception& e) {
    std::cerr << "Exception: " << e.what() << "\n";
    return 1;
  }

  Arrangement arr;
  if (! load_polygons_from_file(filename, arr)) return 1;

  CGAL::Graphics_scene_options<Arrangement, Vertex_const_handle, Halfedge_const_handle, Face_const_handle> gso;
  gso.ignore_all_vertices(true);
  gso.ignore_all_edges(true);
  gso.colored_face = [](const Arrangement&, Face_const_handle) -> bool { return true; };
  gso.face_color = [] (const Arrangement&, Face_const_handle fh) -> CGAL::IO::Color {
    if (fh->is_unbounded()) return CGAL::IO::Color(100, 125, 200);
    return get_random_color(CGAL::get_default_random());
  };
  CGAL::draw(arr, gso, "2D Arrangement");

  // Print results
  std::size_t i = 0;
  for (auto it = arr.faces_begin(); it != arr.faces_end(); ++it) {
    std::cout << "Face (" << i++ << "): ";
    for (const auto& name : it->data()) std::cout << name << " ";
    std::cout << std::endl;
  }
  return 0;
}
