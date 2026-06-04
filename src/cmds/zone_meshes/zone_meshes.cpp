#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <algorithm>
#include <variant>

#include <boost/program_options.hpp>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_polyhedron.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/generators.h>

#include "cgalex/Paths.h"

namespace po = boost::program_options;
namespace fs = std::filesystem;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

#define USE_SURFACE_MESH

using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Plane = Kernel::Plane_3;
using Vector = Kernel::Vector_3;
using Segment = Kernel::Segment_3;
using Line = Kernel::Line_3;
using Ray = Kernel::Ray_3;
#if defined(USE_SURFACE_MESH)
using Mesh = CGAL::Surface_mesh<Point>;
#else
using Mesh = CGAL::Polyhedron_3<Kernel>;
#endif
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits_3<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;
using Segment_intersection = std::optional<Tree::Intersection_and_primitive_id<Segment>::Type>;
using Plane_intersection = std::optional<Tree::Intersection_and_primitive_id<Plane>::Type>;
using Primitive_id = Tree::Primitive_id;

using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

template <typename InputIterator>
void sort_intersections(InputIterator begin, InputIterator end, const Line& line) {
  // Line direction and reference point
  Vector dir = line.to_vector();
  Point ref = line.point();

  // Sort by signed parameter along the line
  std::sort(begin, end,
            [&](Segment_intersection a, Segment_intersection b) {
              const Point* pa = std::get_if<Point>(&(a->first));
              if (! pa) {
                Segment* seg = std::get_if<Segment>(&(a->first));
                CGAL_assertion(seg);
                pa = &(seg->source());
              }
              const Point* pb = std::get_if<Point>(&(b->first));
              if (! pb) {
                Segment* seg = std::get_if<Segment>(&(b->first));
                CGAL_assertion(seg);
                pb = &(seg->source());
              }
              auto ta = ((*pa - ref) * dir);
              auto tb = ((*pb - ref) * dir);
              return ta < tb;
            });
}

template <typename InputIterator>
void print_intersections(InputIterator begin, InputIterator end) {
  for (auto it = begin; it != end; ++it) {
    Segment_intersection x = *it;
    const Point* p = std::get_if<Point>(&(x->first));
    if (p) {
      std::cout << "intersection point " << *p << " with " << x->second << std::endl;
      continue;
    }
    const Segment* s = std::get_if<Segment>(&(x->first));
    CGAL_assertion(s);
    std::cout << "intersection segment " << *s << " with " << x->second << std::endl;
  }
}

//!
template <typename InputIterator, typename NormalMap, typename OutputIterator>
OutputIterator process_intersections(InputIterator begin, InputIterator end, const Line& line,
                                     NormalMap& normals, OutputIterator oi) {
  Vector dir = line.to_vector();
  Point start;
  bool penetrated = false;
  for (auto it = begin; it != end; ++it) {
    Segment_intersection x = *it;

    const Segment* s = std::get_if<Segment>(&(x->first));
    if (s) {
      *oi++ = *s;
      continue;
    }

    const Point* p = std::get_if<Point>(&(x->first));
    CGAL_assertion(p);
    auto next = it;
    ++next;

    if (penetrated) {
      *oi++ = Segment(start, *p);
      penetrated = false;

      // Deduplicate points
      while (next != end) {
        Segment_intersection next_x = *next;
        if (std::get_if<Segment>(&(next_x->first))) break;
        const Point* next_p = std::get_if<Point>(&(next_x->first));
        if (*p != *next_p) break;
        it = next++;
      }
      continue;
    }

    if (next == end) {
      // This must be a tangent point
      *oi++ = *p;
      continue;
    }

    Segment_intersection next_x = *next;
    if (std::get_if<Segment>(&(next_x->first))) {
      // This must be a tangent point
      *oi++ = *p;
    }

    const Point* next_p = std::get_if<Point>(&(next_x->first));
    CGAL_assertion(next_p);
    if (*p != *next_p) {
      // This must be a penetration point
      start = *p;
      penetrated = true;
      continue;
    }

    // Traverse all incident faces...
    auto fd = x->second;
    const auto& normal = get(normals, fd);
    auto d_sign = CGAL::sign(normal * dir);
    bool tangent = false;

    while (next != end) {
      Segment_intersection next_x = *next;
      if (std::get_if<Segment>(&(next_x->first))) break;
      const Point* next_p = std::get_if<Point>(&(next_x->first));
      if (*p != *next_p) break;
      if (! tangent) {
        auto fd = next_x->second;
        const auto& normal = get(normals, fd);
        auto next_d = normal * dir;
        if (d_sign != CGAL::sign(next_d)) tangent = true;
      }
      it = next++;
    }
    if (! tangent) {
      CGAL_assertion(d_sign == CGAL::POSITIVE);
      penetrated = true;
      start = *p;
      continue;
    }
    *oi++ = *p;
  }
  return oi;
}

//! \brief obtains the default value of the input path
const Path def_input_path() {
  static const Path s_def_input_path(".");
  return s_def_input_path;
}

//!
std::vector<fs::path> find_off_files(const fs::path& dir_path) {
  std::vector<fs::path> files;

  try {
    // Validate the path object
    if (!fs::exists(dir_path) || !fs::is_directory(dir_path)) {
      std::cerr << "Error: Provided path does not exist or is not a directory.\n";
      return files;
    }

    // Iterate through the directory using the path object
    for (const auto& entry : fs::directory_iterator(dir_path)) {
      if (entry.is_regular_file()) {
        const auto& ext = entry.path().extension();

        // Directly compare the extension
        if (ext == ".off" || ext == ".OFF") {
          files.push_back(entry.path());
        }
      }
    }
  }
  catch (const fs::filesystem_error& e) {
    std::cerr << "Filesystem error: " << e.what() << '\n';
  }

  return files;
}

//!
int main(int argc, char* argv[]) {
  Path input_path = def_input_path();
  std::size_t verbose_level = 0;

  try {
    // Define options
    po::options_description desc("Allowed options");
    desc.add_options()
      ("input-path,i", po::value<Path>()->default_value({def_input_path()}), "input path")
      ("verbose,v", po::value<std::size_t>(&verbose_level)->implicit_value(0),
       "set verbosity level (0 = quiet)")
      ;

    // This is a place holder for a positional option
    po::positional_options_description p;

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

    input_path = vm["input-path"].as<Path>();
  }
  catch (std::exception& e) {
    std::cout << "Exception: " << e.what() << "\n";
    return 1;
  }

  // constructs segment query
  Point a(-0.2, 0.2, -0.2);
  Point b(1.3, 0.2, 1.3);
  Line line_query(a, b);

  std::cout << "Searching in: " << input_path << "\n";
  std::vector<fs::path> files = find_off_files(input_path);

  Kernel kernel;
  std::vector<Mesh> meshes;
  if (! files.empty()) {
    meshes.reserve(files.size());
    for (const auto& file : files) {
      Mesh mesh;
      if (! PMP::IO::read_polygon_mesh(file.filename().string(), mesh, params::verbose(true).repair_polygon_soup(true))) {
        std::cerr << "Error: could not read mesh from '" << file.filename().string() << "'.\n";
        return 1;
      }
      meshes.push_back(mesh);
    }
  }
  else {
    meshes.reserve(1);
    Point p(1.0, 0.0, 0.0);
    Point q(0.0, 1.0, 0.0);
    Point r(0.0, 0.0, 1.0);
    Point s(0.0, 0.0, 0.0);
    Mesh mesh;

#if defined(USE_SURFACE_MESH)
    CGAL::make_tetrahedron(p, q, r, s, mesh);
#else
    mesh.make_tetrahedron(p, q, r, s);
#endif
    meshes.push_back(mesh);
  }

  // Compute normals
  auto np = params::geom_traits(kernel);
  for (auto& mesh : meshes) {
    // Compute face normals
    auto normals = mesh.add_property_map<face_descriptor, Vector>("f:normals", CGAL::NULL_VECTOR).first;
    PMP::compute_face_normals(mesh, normals, np);
  }

  // Compute intersection
  for (const auto& mesh : meshes) {
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
    CGAL::add_to_graphics_scene(mesh, scene, gso);
    scene.add_segment(a, b, CGAL::Color(255, 0, 0));
    CGAL::draw_graphics_scene(scene);

    // constructs AABB tree
    Tree tree(faces(mesh).first, faces(mesh).second, mesh);

    // computes all intersections with segment query (as pairs object - primitive_id)
    std::vector<Segment_intersection> face_intersections;
    tree.all_intersections(line_query, std::back_inserter(face_intersections));
    sort_intersections(face_intersections.begin(), face_intersections.end(), line_query);
    print_intersections(face_intersections.begin(), face_intersections.end());

    auto normals_opt = mesh.property_map<face_descriptor, Vector>("f:normals");
    CGAL_assertion(normals_opt.has_value());
    auto normals = normals_opt.value();
    std::vector<std::variant<Point, Segment>> intersections;
    process_intersections(face_intersections.begin(), face_intersections.end(), line_query,
                          normals, std::back_inserter(intersections));

    std::cout << "# intersections: " << intersections.size() << std::endl;
    for (const auto& x : intersections) {
      const Point* p = std::get_if<Point>(&x);
      if (p) {
        std::cout << "Tangent point " << *p << std::endl;
        continue;
      }
      const Segment* seg = std::get_if<Segment>(&x);
      CGAL_assertion(seg);
      std::cout << "Intersection segment " << *seg << std::endl;
    }
  }

  return EXIT_SUCCESS;
}
