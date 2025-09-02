#include <iostream>
#include <string>
#include <filesystem>

#include <boost/program_options.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

namespace po = boost::program_options;
namespace fs = std::filesystem;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

//! Not in used yet
template <typename Mesh>
bool is_small_hole(typename boost::graph_traits<Mesh>::halfedge_descriptor h, Mesh& mesh,
                   double max_hole_diam, std::size_t max_num_hole_edges) {
  using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  std::size_t num_hole_edges = 0;
  CGAL::Bbox_3 hole_bbox;
  for (halfedge_descriptor hc : CGAL::halfedges_around_face(h, mesh)) {
    const auto& p = mesh.point(target(hc, mesh));

    hole_bbox += p.bbox();
    ++num_hole_edges;

    // Exit early, to avoid unnecessary traversal of large holes
    if (num_hole_edges > max_num_hole_edges) return false;
    if (hole_bbox.xmax() - hole_bbox.xmin() > max_hole_diam) return false;
    if (hole_bbox.ymax() - hole_bbox.ymin() > max_hole_diam) return false;
    if (hole_bbox.zmax() - hole_bbox.zmin() > max_hole_diam) return false;
  }

  return true;
}

//!
template <typename Mesh>
bool fill_mesh(Mesh& mesh, std::size_t verbose = 0) {
  using vertex_descriptor = typename boost::graph_traits<Mesh>::vertex_descriptor;
  using halfedge_descriptor = typename boost::graph_traits<Mesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<Mesh>::face_descriptor;

  // collect one halfedge per boundary cycle
  std::vector<halfedge_descriptor> border_cycles;
  PMP::extract_boundary_cycles(mesh, std::back_inserter(border_cycles));

  std::size_t nb_holes = 0;
  for (halfedge_descriptor h : border_cycles) {
    // if (max_hole_diam > 0 && max_num_hole_edges > 0 && ! is_small_hole(h, mesh, max_hole_diam, max_num_hole_edges))
    //   continue;

    std::vector<face_descriptor>  patch_facets;
    std::vector<vertex_descriptor> patch_vertices;
    bool success =
      std::get<0>(PMP::triangulate_refine_and_fair_hole(mesh,
                                                        h,
                                                        params::face_output_iterator(std::back_inserter(patch_facets))
                                                        .vertex_output_iterator(std::back_inserter(patch_vertices))));

    if (verbose >= 2) {
      std::cout << "  Number of facets in constructed patch: " << patch_facets.size() << std::endl;
      std::cout << "  Number of vertices in constructed patch: " << patch_vertices.size() << std::endl;
      std::cout << "  Is fairing successful: " << success << std::endl;
    }
    if (! success) return false;
    ++nb_holes;
  }

  if (verbose >= 1) {
    std::cout << nb_holes << " holes have been filled" << std::endl;
  }

  return true;
}

//!
bool check_mesh(const std::string& filename, std::size_t verbose = 0, bool do_draw = false) {
  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

  using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
  using Face_color_map = Mesh::Property_map<Mesh::Face_index, CGAL::IO::Color>;

  using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
  using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
  using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

  Mesh mesh;

  bool has_colors = false;
  Face_color_map fcolors;

  if (! CGAL::IO::read_polygon_mesh(filename, mesh, params::verbose(true).repair_polygon_soup(true))) {
    std::cerr << "Error: could not read mesh from '" << filename << "'.\n";
    return 1;
  }

  auto face_property_maps = mesh.properties<face_descriptor>();
  if (verbose >= 2) {
    for (const auto& pm : face_property_maps) std::cout << pm << std::endl;
    std::cout << filename << " loaded with " << num_vertices(mesh) << " vertices, " << num_faces(mesh) << " faces\n";
  }

  if (do_draw) {
    CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso;
    gso.ignore_all_vertices(true);
    gso.ignore_all_edges(true);
    auto fcolors_optional = mesh.property_map<face_descriptor, CGAL::IO::Color>("f:color");
    if (fcolors_optional) {
      has_colors = true;
      fcolors = fcolors_optional.value();
    }
    if (has_colors) {
      gso.colored_face=[](const Mesh&, face_descriptor) -> bool { return true; };
      gso.face_color=[&](const Mesh& mesh, face_descriptor fd) -> CGAL::IO::Color { return get(fcolors, fd); };
    }
    CGAL::draw(mesh, gso, filename.c_str());
  }

  // Closed
#if 0
  auto is_closed = CGAL::is_closed(mesh);
  if (! is_closed) {
    std::cerr << "Warning: " << filename << " is not closed; filling...";
    auto success = fill_mesh(mesh);
    if (success) {
      std::cout << "done\n";
      is_closed = true;
      // CGAL::IO::write_polygon_mesh("out.off", mesh, params::stream_precision(17));
      // std::cout << "Mesh written to: out.off" << std::endl;
    }
    else {
      std::cerr << "failed\n";
      return false;
    }
  }

  // Now try again
#endif
  auto is_closed = CGAL::is_closed(mesh);
  if (! is_closed) {
    std::cerr << filename << " is not closed\n";
    return false;
  }

  // Triangular mesh
  auto is_tri = CGAL::is_triangle_mesh(mesh);
  if (! is_tri) {
    std::cerr << filename << " is not triangular\n";
    return false;
  }

  // Connected components
  auto fccmap = mesh.add_property_map<face_descriptor, std::size_t>("f:CC").first;
  std::size_t num_ccs = PMP::connected_components(mesh, fccmap);
  if (num_ccs > 1) {
    std::cout << filename << " number of connected components: " << num_ccs << "\n";
  }

  // Does bound a volume
  CGAL_assertion(is_closed && is_tri);
  std::vector<bool> isoo;
  bool does_bound_volume =
    PMP::does_bound_a_volume(mesh, params::is_cc_outward_oriented(std::reference_wrapper(isoo)));
  if (! does_bound_volume) {
    std::cerr << filename << " does not bound a volume\n";
    return false;
  }
  for (auto ccid = 0; ccid < num_ccs; ++ccid) {
    if (! isoo[ccid]) {
      std::cerr << filename << " connected component " << ccid << " is not outward oriented\n";
      return false;
    }
  }
  return true;
}

//!
int main(int argc, char* argv[]) {
  try {
    // double max_hole_diam = (argc > 2) ? boost::lexical_cast<double>(argv[2]): -1.0;
    // std::size_t max_num_hole_edges = (argc > 3) ? boost::lexical_cast<int>(argv[3]) : -1;
    bool do_draw = false;
    std::size_t verbose = 0;
    std::string input_dir = ".";

    // Define options
    po::options_description desc("Allowed options");
    desc.add_options()
      ("help,h", "produce help message")
      ("verbose,v", po::value<std::size_t>(&verbose)->implicit_value(1),
       "set verbosity level (0 = quiet, 1 = show directory, 2 = also count files)")
      ("draw,d", po::value<bool>(&do_draw)->implicit_value(false), "set draw")
      ("input-directory,i", po::value<std::string>(&input_dir)->default_value("."),
       "input directory (default: current directory)");

    // Parse options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
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

    // Traverse .off files
    std::size_t count = 0;
    if (fs::exists(input_dir) && fs::is_directory(input_dir)) {
      for (const auto& entry : fs::directory_iterator(input_dir)) {
        if (entry.is_regular_file() && entry.path().extension() == ".off") {
          if (verbose >= 1) std::cout << entry.path().string() << ": checking\n";
          check_mesh(entry.path().string(), verbose, do_draw);
          ++count;
        }
      }
    }
    else {
      std::cerr << "Error: '" << input_dir << "' is not a valid directory.\n";
      return 1;
    }

    // Verbosity: show number of files found
    if (verbose >= 1) {
      std::cout << "Found " << count << " .off file(s)\n";
    }
  }
  catch (std::exception& e) {
    std::cerr << "Exception: " << e.what() << "\n";
    return 1;
  }

  return 0;
}
