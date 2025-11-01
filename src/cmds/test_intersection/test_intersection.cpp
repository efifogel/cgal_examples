#include <iostream>
#include <list>

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
using Ray = Kernel::Ray_3;
#if defined(USE_SURFACE_MESH)
using Mesh = CGAL::Surface_mesh<Point>;
#else
using Mesh = CGAL::Polyhedron_3<Kernel>;
#endif
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Mesh>;
using Traits = CGAL::AABB_traits_3<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;
using Segment_intersection = std::optional< Tree::Intersection_and_primitive_id<Segment>::Type>;
using Plane_intersection = std::optional< Tree::Intersection_and_primitive_id<Plane>::Type>;
using Primitive_id = Tree::Primitive_id;

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
  Segment segment_query(a, b);

  // tests intersections with segment query
  if (tree.do_intersect(segment_query)) std::cout << "intersection(s)" << std::endl;
  else std::cout << "no intersection" << std::endl;

  // computes #intersections with segment query
  std::cout << tree.number_of_intersected_primitives(segment_query)
            << " intersection(s)" << std::endl;

  // computes first encountered intersection with segment query
  // (generally a point)
  Segment_intersection intersection = tree.any_intersection(segment_query);
  if (intersection) {
    // gets intersection object
    Point* p = std::get_if<Point>(&(intersection->first));
    if (p) std::cout << "intersection object is a point " << *p << "with face "<< intersection->second  <<  std::endl;
  }

  // computes all intersections with segment query (as pairs object - primitive_id)
  std::list<Segment_intersection> intersections;
  tree.all_intersections(segment_query, std::back_inserter(intersections));
  for (auto x : intersections) {
    Point* p = std::get_if<Point>(&(x->first));
    if (p) {
      std::cout << "intersection point " << *p << " with " << x->second << std::endl;
      continue;
    }
    Segment* s = std::get_if<Segment>(&(x->first));
    if (s) {
      std::cout << "intersection segment " << *s << " with " << x->second << std::endl;
      continue;
    }
  }

  // computes all intersected primitives with segment query as primitive ids
  std::list<Primitive_id> primitives;
  tree.all_intersected_primitives(segment_query, std::back_inserter(primitives));

  // constructs plane query
  Point base(0.0, 0.0, 0.5);
  Vector vec(0.0, 0.0, 1.0);
  Plane plane_query(base,vec);

  // computes first encountered intersection with plane query
  // (generally a segment)
  Plane_intersection plane_intersection = tree.any_intersection(plane_query);
  if (plane_intersection){
    if (std::get_if<Segment>(&(plane_intersection->first))){
      Segment* s = std::get_if<Segment>(&(plane_intersection->first));
      std::cout << "one intersection object is the segment " << *s << "with face "<< intersection->second  <<  std::endl;
    }
  }

  return EXIT_SUCCESS;
}
