#include <iostream>
#include <vector>
#include <variant>

#define USE_SURFACE_MESH

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/draw_polyhedron.h>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/generators.h>

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
void sort_intersections(const Line& line, InputIterator begin, InputIterator end) {
  // Line direction and reference point
  Vector dir = line.to_vector();
  Point ref = line.point();

  // Sort by signed parameter along the line
  std::sort(begin, end,
            [&](Segment_intersection a, Segment_intersection b) {
              const Point* pa = std::get_if<Point>(&(a->first));
              if (not pa) {
                Segment* seg = std::get_if<Segment>(&(a->first));
                CGAL_assertion(seg);
                pa = &(seg->source());
              }
              const Point* pb = std::get_if<Point>(&(b->first));
              if (not pb) {
                Segment* seg = std::get_if<Segment>(&(b->first));
                CGAL_assertion(seg);
                pb = &(seg->source());
              }
              auto ta = ((*pa - ref) * dir);
              auto tb = ((*pb - ref) * dir);
              return ta > tb;
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

template <typename InputIterator, typename OutputIterator>
OutputIterator process_intersections(InputIterator begin, InputIterator end, OutputIterator oi) {
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
    std::cout << "XXXX current point " << *p << std::endl;
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
    std::cout << "XXXX next point " << *next_p << std::endl;
    CGAL_assertion(next_p);
    if (*p != *next_p) {
      // This must be a penetration point
      start = *p;
      penetrated = true;
      continue;
    }

    // Traverse all incident faces...
  }
  return oi;
}

int main(int argc, char* argv[]) {
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

  // constructs AABB tree
  Tree tree(faces(mesh).first, faces(mesh).second, mesh);

  // constructs segment query
  Point a(-0.2, 0.2, -0.2);
  Point b(1.3, 0.2, 1.3);
  Line line_query(a, b);

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

  // computes all intersections with segment query (as pairs object - primitive_id)
  std::vector<Segment_intersection> intersections;
  tree.all_intersections(line_query, std::back_inserter(intersections));
  sort_intersections(line_query, intersections.begin(), intersections.end());
  print_intersections(intersections.begin(), intersections.end());

  std::vector<std::variant<Point, Segment>> result;
  process_intersections(intersections.begin(), intersections.end(), std::back_inserter(result));
  std::cout << "# intersections: " << result.size() << std::endl;
  for (const auto& x : result) {
    const Point* p = std::get_if<Point>(&x);
    if (p) {
      std::cout << "Tangent point " << *p << std::endl;
      continue;
    }
    const Segment* seg = std::get_if<Segment>(&x);
    CGAL_assertion(seg);
    std::cout << "Intersection segment " << *seg << std::endl;
  }

  return EXIT_SUCCESS;
}
