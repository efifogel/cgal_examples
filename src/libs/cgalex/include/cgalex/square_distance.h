#ifndef CGALEX_SQUARE_DISTANCE_H
#define CGALEX_SQUARE_DISTANCE_H

#include <CGAL/enum.h>

/*! Given a halfedge circulator around a vertex of a Gaussian map, compute the
 * sqaure distance between the origin and the facet of the primal polyhedron
 * that is dual to the vertex.
 */
template <typename HalfedgeAroundVertexCirculator, typename Kernel>
typename Kernel::FT square_distance(HalfedgeAroundVertexCirculator hec, const Kernel& kernel) {
  const auto& p1 = hec++->face()->data();
  const auto& p2 = hec++->face()->data();
  const auto& p3 = hec->face()->data();
  auto plane = kernel.construct_plane_3_object()(p1, p2, p3);
  typename Kernel::Point_3 origin(CGAL::ORIGIN);
  return CGAL::squared_distance(origin, plane);
}

/*! Given a halfedge of a Gaussian map and a direction, compute the sqaure
 * distance between the origin and the point where the line in the given
 * direction emenating from the origin intersects the facet of the primal
 * polyhedron that is dual to the halfedge endpoint.
 */
template <typename HalfedgeHandle, typename Kernel>
typename Kernel::FT
square_distance(const typename Kernel::Direction_3& direction, HalfedgeHandle he, const Kernel& kernel) {
  const auto& q = he->face()->data();
  auto plane = kernel.construct_plane_3_object()(q, direction);
  typename Kernel::Point_3 origin(CGAL::ORIGIN);
  return CGAL::squared_distance(origin, plane);
}

#endif
