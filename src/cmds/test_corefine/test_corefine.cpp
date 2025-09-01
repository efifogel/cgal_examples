#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

// #include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Basic_viewer.h>
#include <CGAL/draw_surface_mesh.h>

// using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
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

  Mesh mesh1, mesh2;
  if (! PMP::IO::read_polygon_mesh(filename1, mesh1) || ! PMP::IO::read_polygon_mesh(filename2, mesh2)) {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }

  CGAL::Graphics_scene scene;
  CGAL::add_to_graphics_scene(mesh1, scene);
  CGAL::add_to_graphics_scene(mesh2, scene);
  CGAL::draw_graphics_scene(scene);

  // create a property on edges to indicate whether they are constrained
  Mesh::Property_map<edge_descriptor,bool> is_constrained_map =
    mesh1.add_property_map<edge_descriptor,bool>("e:is_constrained", false).first;

  // update mesh1 to contain the mesh bounding the difference
  // of the two input volumes.
  bool valid_difference =
    PMP::corefine_and_compute_difference(mesh1,
                                         mesh2,
                                         mesh1,
                                         params::default_values(), // default parameters for mesh1
                                         params::default_values(), // default parameters for mesh2
                                         params::edge_is_constrained_map(is_constrained_map));

  if (! valid_difference) {
    std::cerr << "Difference could not be computed\n";
    return 1;
  }

  std::cout << "Difference was successfully computed\n";
  CGAL::draw(mesh1, "difference");
  CGAL::IO::write_polygon_mesh("difference.off", mesh1, CGAL::parameters::stream_precision(17));

#if 0
  // collect faces incident to a constrained edge
  std::vector<face_descriptor> selected_faces;
  std::vector<bool> is_selected(num_faces(mesh1), false);
  for (edge_descriptor e : edges(mesh1))
    if (is_constrained_map[e]) {
      // insert all faces incident to the target vertex
      for (halfedge_descriptor h : halfedges_around_target(halfedge(e,mesh1),mesh1)) {
        if (! is_border(h, mesh1)) {
          face_descriptor f=face(h, mesh1);
          if (! is_selected[f]) {
            selected_faces.push_back(f);
            is_selected[f]=true;
          }
        }
      }
    }

  // increase the face selection
  CGAL::expand_face_selection(selected_faces, mesh1, 2,
    Vector_pmap_wrapper(is_selected), std::back_inserter(selected_faces));

  std::cout << selected_faces.size()
            << " faces were selected for the remeshing step\n";

  // remesh the region around the intersection polylines
  PMP::isotropic_remeshing(selected_faces, 0.02, mesh1,
                           params::edge_is_constrained_map(is_constrained_map));

  CGAL::IO::write_polygon_mesh("difference_remeshed.off", mesh1, CGAL::parameters::stream_precision(17));
#endif

  return 0;
}
