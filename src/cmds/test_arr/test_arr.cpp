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
#include <CGAL/Arr_observer.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/draw_arrangement_2.h>
#include <CGAL/Graphics_scene_options.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Traits = CGAL::Arr_segment_traits_2<Kernel>;
using Name = std::string;
using Names = std::list<std::string>;
using Dcel = CGAL::Arr_extended_dcel<Traits, int, Names, std::pair<bool, Names>>;
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

// An observer that receives notifications of face splits and face mergers.
class My_observer : public CGAL::Arr_observer<Arrangement> {
private:
  std::string m_name;

public:
  My_observer(const std::string& name) : m_name(name) {}

  My_observer(Arrangement& arr) : CGAL::Arr_observer<Arrangement>(arr) {}
  virtual void after_create_edge(Halfedge_handle e) {
    e->data().push_back(m_name);
    e->twin()->data().push_back(m_name);
  }

  virtual void 	after_split_edge(Halfedge_handle e1, Halfedge_handle e2) {
    e2->set_data(e1->data());
    e2->twin()->set_data(e1->data());
  }
};

//! Insert a polygon into the arrangement
template <typename InputIterator>
void add_polygon(Arrangement& arr, const std::string& name, InputIterator begin, InputIterator end) {
  std::cout << "Adding polygon " << name << "\n";
  My_observer observer(name);
  observer.attach(arr);
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
  observer.detach();
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

    add_polygon(arr, name, points.begin(), points.end());
  }

  return true;
}

//!
Names update_names(Halfedge_handle e, Names& names, const Traits& traits) {
  auto cmp_endpoints = traits.compare_endpoints_xy_2_object();
  Names new_names;
  auto& delta_names = e->data();
  auto he_dir = e->direction();
  auto curve_dir = cmp_endpoints(e->curve());
  if (((curve_dir == CGAL::LARGER) && (he_dir == CGAL::ARR_LEFT_TO_RIGHT)) ||
      ((curve_dir == CGAL::SMALLER) && (he_dir == CGAL::ARR_RIGHT_TO_LEFT))) {
    new_names = names;
    new_names.insert(new_names.end(), delta_names.begin(), delta_names.end());
    return new_names;
  }
  for (const auto& name : names) {
    if (std::find(delta_names.begin(), delta_names.end(), name) == delta_names.end()) new_names.push_back(name);
  }
  return new_names;
}

//!
void process(Face_handle f, Names& names, const Traits& traits) {
  if (f->data().first) return;
  f->set_data(std::make_pair(true, names));

  for (auto it = f->inner_ccbs_begin(); it != f->inner_ccbs_end(); ++it) {
    auto curr = *it;
    do {
      auto inner_face = curr->twin()->face();
      if (inner_face->data().first) continue;
      Names new_names = update_names(curr, names, traits);
      process(inner_face, new_names, traits);
    } while (++curr != *it);
  }

  // Traverse the outer boundary.
  if (f->is_unbounded()) return;

  for (auto it = f->outer_ccbs_begin(); it != f->outer_ccbs_end(); ++it) {
    auto curr = *it;
    do {
      auto outer_face = curr->twin()->face();
      if (outer_face->data().first) continue;
      Names new_names = update_names(curr, names, traits);
      process(outer_face, new_names, traits);
    } while (++curr != *it);
  }
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

  // Initialize the "visited" flag of all faces
  for (auto it = arr.faces_begin(); it != arr.faces_end(); ++it) it->data().first = false;

  // Traverse the arrangement and for each face assign the names of the covering polygons
  const Traits& traits = *(arr.geometry_traits());
  Names names;
  process(arr.unbounded_face(), names, traits);

  // Print results
  std::size_t i = 0;
  for (auto it = arr.faces_begin(); it != arr.faces_end(); ++it) {
    std::cout << "Face (" << i++ << "): ";
    for (const auto& name : it->data().second) std::cout << name << " ";
    std::cout << std::endl;
  }
  return 0;
}
