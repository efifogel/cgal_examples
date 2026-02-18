#ifndef EXTENDED_POLYHEDRON_3_HPP
#define EXTENDED_POLYHEDRON_3_HPP

#include <CGAL/Polyhedron_items_3.h>

template <typename Refs, typename Point>
class Extended_vertex :
  public CGAL::HalfedgeDS_vertex_base<Refs, CGAL::Tag_true, Point> {
private:
  using True = CGAL::Tag_true;
  using Base = typename CGAL::HalfedgeDS_vertex_base<Refs, True, Point>;
  bool m_flag;

public:
  Extended_vertex() : m_flag(false) {}
  Extended_vertex(Point const& p) : Base(p), m_flag(false) {}
  Extended_vertex(Point const& p, bool flag) : Base(p), m_flag(flag) {}
  void set_flag(bool flag) { m_flag = flag; }
  bool flag() const { return m_flag; }
};

struct Extended_polyhedron_items : public CGAL::Polyhedron_items_3 {
  using True = CGAL::Tag_true;
  template <typename Refs, typename Traits>
  struct Vertex_wrapper
  { using Vertex = Extended_vertex<Refs, typename Traits::Point_3>; };

  template <typename Refs, typename Traits>
  struct Halfedge_wrapper
  { using Halfedge = CGAL::HalfedgeDS_halfedge_base<Refs, True,True,True>; };

  template <typename Refs, typename Traits>
  struct Face_wrapper {
    using Face =
      CGAL::HalfedgeDS_face_base<Refs, True, typename Traits::Plane_3>;
  };
};

#endif
