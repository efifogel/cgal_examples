#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <algorithm>
#include <variant>
#include <ranges>

#include <boost/program_options.hpp>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Arrangement_on_curve_1/Arrangement_on_curve_1.h>
#include <CGAL/Arrangement_on_curve_1/Unbounded_topology_traits.h>
#include <CGAL/Arrangement_on_curve_1/overlay.h>
#include <CGAL/Arrangement_on_curve_1/Line_3_traits_1.h>
#include <CGAL/Arrangement_on_curve_1/insert.h>
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
namespace Aoc1 = CGAL::Arrangement_on_curve_1;

#define USE_SURFACE_MESH

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
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

using Geometry_traits_1 = Aoc1::Line_3_traits_1<Kernel>;
using Point_1 = Geometry_traits_1::Point_1;
using Topology_traits = Aoc1::Unbounded_topology_traits<Point_1, std::vector<std::size_t>, std::vector<std::size_t>>;
using Arrangement1 = Aoc1::Arrangement_on_curve_1<Geometry_traits_1, Topology_traits>;

template <typename ArrangementA, typename ArrangementB, typename ArrangementR>
class Combine_overlay_observer {
public:
  using Vertex_const_descriptor_a = typename ArrangementA::Vertex_const_descriptor;
  using Vertex_const_descriptor_b = typename ArrangementB::Vertex_const_descriptor;
  using Vertex_descriptor_r = typename ArrangementR::Vertex_descriptor;

  using Edge_const_descriptor_a = typename ArrangementA::Edge_const_descriptor;
  using Edge_const_descriptor_b = typename ArrangementB::Edge_const_descriptor;
  using Edge_descriptor_r = typename ArrangementR::Edge_descriptor;

  using Vertex_data_map_a = typename ArrangementA::Topology_traits::Vertex_data_map;
  using Vertex_data_map_b = typename ArrangementB::Topology_traits::Vertex_data_map;
  using Vertex_data_map_r = typename ArrangementR::Topology_traits::Vertex_data_map;

  using Edge_data_map_a = typename ArrangementA::Topology_traits::Edge_data_map;
  using Edge_data_map_b = typename ArrangementB::Topology_traits::Edge_data_map;
  using Edge_data_map_r = typename ArrangementR::Topology_traits::Edge_data_map;

private:
  Vertex_data_map_r m_v_map_r;
  Edge_data_map_r m_e_map_r;

  Vertex_data_map_a m_v_map_a;
  Edge_data_map_a m_e_map_a;

  Vertex_data_map_b m_v_map_b;
  Edge_data_map_b m_e_map_b;

  // Helper utility to combine two vectors of indices
  std::vector<std::size_t> combine_vectors(const std::vector<std::size_t>& vec_a,
                                           const std::vector<std::size_t>& vec_b) const {
    std::vector<std::size_t> combined;
    combined.reserve(vec_a.size() + vec_b.size());
    combined.insert(combined.end(), vec_a.begin(), vec_a.end());
    combined.insert(combined.end(), vec_b.begin(), vec_b.end());
    return combined;
  }

public:
  Combine_overlay_observer(const ArrangementA& arr_a, const ArrangementB& arr_b, ArrangementR& arr_r) :
    m_v_map_r(arr_r.vertex_data_map()),
    m_e_map_r(arr_r.edge_data_map()),
    m_v_map_a(arr_a.vertex_data_map()),
    m_e_map_a(arr_a.edge_data_map()),
    m_v_map_b(arr_b.vertex_data_map()),
    m_e_map_b(arr_b.edge_data_map())
  {}

  // ============================================================================
  // 1. Two vertices coincide
  // ============================================================================
  void create_vertex(Vertex_const_descriptor_a v_a, Vertex_const_descriptor_b v_b, Vertex_descriptor_r v_res) {
    const auto& vec_a = get(m_v_map_a, v_a);
    const auto& vec_b = get(m_v_map_b, v_b);
    put(m_v_map_r, v_res, combine_vectors(vec_a, vec_b));
  }

  // ============================================================================
  // 2. Vertex A splits Edge B
  // ============================================================================
  void create_vertex(Vertex_const_descriptor_a v_a, Edge_const_descriptor_b e_b, Vertex_descriptor_r v_res) {
    const auto& vec_a = get(m_v_map_a, v_a);
    const auto& vec_b = get(m_e_map_b, e_b);
    put(m_v_map_r, v_res, combine_vectors(vec_a, vec_b));
  }

  // ============================================================================
  // 3. Edge A is split by Vertex B
  // ============================================================================
  void create_vertex(Edge_const_descriptor_a e_a, Vertex_const_descriptor_b v_b, Vertex_descriptor_r v_res) {
    const auto& vec_a = get(m_e_map_a, e_a);
    const auto& vec_b = get(m_v_map_b, v_b);
    put(m_v_map_r, v_res, combine_vectors(vec_a, vec_b));
  }

  // ============================================================================
  // 4. Two edges overlap over a shared interval
  // ============================================================================
  void create_edge(Edge_const_descriptor_a e_a, Edge_const_descriptor_b e_b, Edge_descriptor_r e_res) {
    const auto& vec_a = get(m_e_map_a, e_a);
    const auto& vec_b = get(m_e_map_b, e_b);
    put(m_e_map_r, e_res, combine_vectors(vec_a, vec_b));
  }
};

template <typename InputIterator>
void sort_intersections(InputIterator begin, InputIterator end, const Line& line) {
  // Line direction and reference point
  Vector dir = line.to_vector();
  Point ref = line.point();

  // Sort by signed parameter along the line
  std::sort(begin, end,
            [&](const Segment_intersection& a, const Segment_intersection& b) {
              const Point* pa_ptr = std::get_if<Point>(&(a->first));

              // This reference extends the lifetime of the temporary if source() is called
              const Point& pa_ref = pa_ptr ? *pa_ptr : std::get_if<Segment>(&(a->first))->source();

              // Extract or construct Point B
              const Point* pb_ptr = std::get_if<Point>(&(b->first));
              const Point& pb_ref = pb_ptr ? *pb_ptr : std::get_if<Segment>(&(b->first))->source();

              // Both references are 100% valid here
              auto ta = ((pa_ref - ref) * dir);
              auto tb = ((pb_ref - ref) * dir);
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
template <typename GeometryTraits, typename TopologyTraits>
void traverse_left_to_right(const Aoc1::Arrangement_on_curve_1<GeometryTraits, TopologyTraits>& arr) {
  // If the arrangement has no vertices, it contains exactly one completely unbounded edge
  if (arr.is_empty()) {
    std::cout << "The arrangement is empty (contains 1 fully unbounded line edge).\n";
    return;
  }

  const auto& vertices_range = arr.vertices();
  auto p_map = arr.vertex_point_map();
  auto v_data_map = arr.vertex_data_map();
  auto e_data_map = arr.edge_data_map();

  std::cout << "Starting Traversal:\n";

  // Start with the leftmost unbounded edge (-inf, ...)
  auto curr_e = arr.unbounded_left_edge();
  std::cout << "  [Unbounded Left Edge: ";
  std::ranges::copy(get(e_data_map, curr_e), std::ostream_iterator<std::size_t>(std::cout, " "));
  std::cout << "]\n";

  // Proceed to the right by traversing the topological graph
  while (arr.has_right_vertex(curr_e)) {
    auto curr_v = arr.right_vertex(curr_e);
    // Process the current vertex
    std::cout << "  -> Vertex: (" << get(p_map, curr_v) << "), ";
    std::ranges::copy(get(v_data_map, curr_v), std::ostream_iterator<std::size_t>(std::cout, " "));
    std::cout << "\n";

    // Move to its right incident edge
    curr_e = arr.right_edge(curr_v);

    // If this edge does not have a target/right vertex, we have reached the most right edge
    if (! arr.has_right_vertex(curr_e)) {
      std::cout << "  -> [Unbounded Right Edge: ";
      std::ranges::copy(get(e_data_map, curr_e), std::ostream_iterator<std::size_t>(std::cout, " "));
      std::cout << "]\n";
      break;
    }
    else {
      std::cout << "  -> Edge: ";
      std::ranges::copy(get(e_data_map, curr_e), std::ostream_iterator<std::size_t>(std::cout, " "));
      std::cout << "\n";
    }
  }

  std::cout << "Traversal complete.\n";
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

  auto traits_ptr = std::make_shared<const Geometry_traits_1>(line_query);
  Arrangement1 arr(traits_ptr);

  // Compute intersection
  std::size_t i = 0;
  for (const auto& mesh : meshes) {
    CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso;
    gso.ignore_all_vertices(true);
    gso.ignore_all_edges(true);
    gso.colored_face = [](const Mesh&, typename boost::graph_traits<Mesh>::face_descriptor) -> bool { return true; };
    gso.face_color = [] (const Mesh&, typename boost::graph_traits<Mesh>::face_descriptor fh) -> CGAL::IO::Color {
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
    // print_intersections(face_intersections.begin(), face_intersections.end());

    auto normals_opt = mesh.property_map<face_descriptor, Vector>("f:normals");
    CGAL_assertion(normals_opt.has_value());
    auto normals = normals_opt.value();
    std::vector<std::variant<Point, Segment>> intersections;
    process_intersections(face_intersections.begin(), face_intersections.end(), line_query,
                          normals, std::back_inserter(intersections));

    Arrangement1 arr_mesh(traits_ptr);

    // Retrieve the read/write property maps from the arrangement facade
    auto v_data_map = arr_mesh.vertex_data_map();
    auto e_data_map = arr_mesh.edge_data_map();

    for (const auto& x : intersections) {
      const Point* p = std::get_if<Point>(&x);
      if (p) {
        // std::cout << "Tangent point " << *p << std::endl;
        auto v = Aoc1::insert(arr_mesh, *p);
        get(v_data_map, v).push_back(i);
        continue;
      }
      const Segment* seg = std::get_if<Segment>(&x);
      CGAL_assertion(seg);
      // std::cout << "Intersection segment " << *seg << std::endl;

      // Fetch the segment endpoints
      Point_1 p_src = seg->source();
      Point_1 p_tgt = seg->target();

      // Ensure we insert them in sorted order along the 1D line track
      auto comp = arr_mesh.geometry_traits_1().compare_x_1_object();
      if (comp(p_src, p_tgt) == CGAL::LARGER) std::swap(p_src, p_tgt);

      // Insert the two endpoints into the 1D arrangement
      // Since no vertices exist in this range, these operations split the sequence cleanly
      auto v_first = CGAL::Arrangement_on_curve_1::insert(arr_mesh, p_src);
      auto v_second = CGAL::Arrangement_on_curve_1::insert(arr_mesh, p_tgt);

      // Append index 'i' to the vertex data vectors
      get(v_data_map, v_first).push_back(i);
      get(v_data_map, v_second).push_back(i);

      // Locate the unique bounded edge connecting v_first and v_second
      // In our 1D topology, the edge to the right of v_first spans precisely to v_second
      auto e_bet = arr_mesh.right_edge(v_first);

      // Append index 'i' to the connecting edge's data vector
      get(e_data_map, e_bet).push_back(i);
    }
    Arrangement1 arr_tmp(traits_ptr);
    using Observer = Combine_overlay_observer<Arrangement1, Arrangement1, Arrangement1>;
    Observer observer(arr, arr_mesh, arr_tmp);
    Aoc1::overlay(arr, arr_mesh, arr_tmp, observer);
    std::swap(arr, arr_tmp);
  }
  traverse_left_to_right(arr);

  return EXIT_SUCCESS;
}
