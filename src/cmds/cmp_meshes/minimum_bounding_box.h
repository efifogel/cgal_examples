#ifndef MINIMUM_BOUNDING_BOX_H
#define MINIMUM_BOUNDING_BOX_H

#include <array>

#include <CGAL/Arr_naive_point_location.h>

#include "cgalex/find_square_width.h"
#include "cgalex/square_distance.h"

template <typename Arrangement_, typename Kernel_>
std::pair<std::array<typename Kernel_::FT, 3>, std::array<typename Kernel_::Direction_3, 3>>
minimum_bounding_box(const Arrangement_& gm, const Kernel_& kernel) {
  using Arrangement = Arrangement_;
  using Kernel = Kernel_;

  using FT = typename Kernel::FT;
  using Point_3 = typename Kernel::Point_3;
  using Direction_3 = typename Kernel::Direction_3;
  using Line_3 = typename Kernel::Line_3;
  using Vector_3 = typename Kernel::Vector_3;

  std::array<FT, 3> mbb_sizes;
  std::array<Direction_3, 3> mbb_directions;

  using Vertex_const_handle = typename Arrangement::Vertex_const_handle;
  using Halfedge_const_handle = typename Arrangement::Halfedge_const_handle;
  using Face_const_handle = typename Arrangement::Face_const_handle;

  auto ctr_plane = kernel.construct_plane_3_object();
  auto ctr_line = kernel.construct_line_3_object();
  auto ctr_cross_product = kernel.construct_cross_product_vector_3_object();
  auto intersect = kernel.intersect_3_object();

  // 1st dimension
  std::tie(mbb_sizes[0], mbb_directions[0]) = find_square_width(gm, kernel);

  // 2nd dimension
  auto plane1 = ctr_plane(CGAL::ORIGIN, mbb_directions[0]);
  mbb_sizes[1] = 0;
  auto cmp1 = [&](const Point_3& p) {
    Vector_3 v(CGAL::ORIGIN, p);
    auto sw = v.squared_length();
    if (sw < mbb_sizes[1]) {
      mbb_sizes[1] = sw;
      mbb_directions[1] = v.direction();
    }
  };
  std::function<void(const Point_3&)> set1;
  auto init1 = [&](const Point_3& p) {
    set1 = cmp1;
    Vector_3 v(CGAL::ORIGIN, p);
    mbb_sizes[1] = v.squared_length();
    mbb_directions[1] = v.direction();
  };
  set1 = init1;

  std::for_each(gm.vertices_begin(), gm.vertices_end(),
                [&](const typename Arrangement::Vertex& vertex) {
                  if (vertex.degree() < 3) return;
                  const auto& dir = vertex.point();
                  auto e(vertex.incident_halfedges());
                  const auto& q = e->face()->data();
                  auto plane = ctr_plane(q, dir);
                  auto res = intersect(plane, plane1);
                  if (! res) return;
                  const auto* lp = std::get_if<Line_3>(&*res);
                  if (! lp) return;
                  auto p = lp->projection(CGAL::ORIGIN);
                  set1(p);
                });

  // 3rd dimension
  mbb_sizes[2] = 0;

  using Naive_pl = CGAL::Arr_naive_point_location<Arrangement>;
  Naive_pl pl(gm);

  auto v0 = mbb_directions[0].vector();
  auto v1 = mbb_directions[1].vector();
  auto v2 = ctr_cross_product(v0, v1);
  mbb_directions[2] = v2.direction();
  auto ctr_point = gm.geometry_traits()->construct_point_2_object();
  auto p = ctr_point(mbb_directions[2]);
  auto obj = pl.locate(p);

  const Face_const_handle* fp;
  if ((fp = std::get_if<Face_const_handle>(&obj))) { // On a face
    auto it = ((*fp)->has_outer_ccb()) ?
      (*fp)->outer_ccbs_begin() : (*fp)->inner_ccbs_begin();
    auto first = *it;
    auto e = first;
    mbb_sizes[2] = square_distance(mbb_directions[2], e, kernel);
    while (++e != first) {
      auto sm = square_distance(mbb_directions[2], e, kernel);
      if (sm < mbb_sizes[2]) mbb_sizes[2] = sm;
    }
    return std::make_pair(mbb_sizes, mbb_directions);
  }

  const Halfedge_const_handle* ep;
  if ((ep = std::get_if<Halfedge_const_handle>(&obj))) { // on an edge
    mbb_sizes[2] = square_distance(mbb_directions[2], (*ep), kernel);
    auto sm = square_distance(mbb_directions[2], (*ep)->twin(), kernel);
    if (sm < mbb_sizes[2]) mbb_sizes[2] = sm;
    return std::make_pair(mbb_sizes, mbb_directions);
  }

  const Vertex_const_handle* vp;
  if ((vp = std::get_if<Vertex_const_handle>(&obj))) { // on a vertex
    // Compute the distance between the origin and the point where the ray in
    // this direction emenating from the origin intersects the primary facet.
    auto e = (*vp)->incident_halfedges();
    mbb_sizes[2] = square_distance(mbb_directions[2], e, kernel);
    return std::make_pair(mbb_sizes, mbb_directions);
  }

  CGAL_error_msg("Invalid object.");
  return std::make_pair(mbb_sizes, mbb_directions);
}

#endif
