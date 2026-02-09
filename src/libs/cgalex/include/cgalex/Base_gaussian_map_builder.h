#ifndef CGALEX_BASE_GAUSSIAN_MAP_BUILDER_H
#define CGALEX_BASE_GAUSSIAN_MAP_BUILDER_H

#include <list>

#include "cgalex/arr_gaussian_map.h"
#include "cgalex/make_x_monotone.h"

class Base_gaussian_map_builder {
public:
  /*! Construct
   */
  Base_gaussian_map_builder(Arrangement& gm) : m_gm(gm) {}

protected:
  // Data memebers
  Arrangement& m_gm;

  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the GM. Each normal defines an end point of the great arc.
   * \param normal1 represents the source normal.
   * \param normal2 represents the target normal.
   * \return the handle for the halfedge directed from the endpoint
   * represented by normal1 toward the endpoint represented by normal2
   * \pre the Gaussian map is empty.
   */
  template <typename Vector_3, typename OutputIterator>
  OutputIterator insert(const Vector_3& normal1, const Vector_3& normal2,
                        OutputIterator oi);

  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the GM. Each normal defines an end point of the great arc.
   * \param normal1 represents the source normal.
   * \param vertex1 the handle of the vertex that is the source of the arc
   * \param normal2 represents the target normal.
   * \return the handle for the halfedge directed from the endpoint
   * represented by normal1 toward the endpoint represented by normal2
   * \param vertex the handle of the vertex that is the source of the arc
   */
  template <typename Vector_3, typename OutputIterator>
  OutputIterator insert(const Vector_3& normal1, Vertex_handle vertex1,
                        const Vector_3& normal2, OutputIterator oi);

  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the GM. Each normal defines an end point of the great arc.
   * \param normal1 represents the source normal.
   * \param normal2 represents the target normal.
   * \param vertex2 the handle of the vertex that is the target of the arc
   * \return the handle for the halfedge directed from the endpoint
   * represented by normal1 toward the endpoint represented by normal2
   */
  template <typename Vector_3, typename OutputIterator>
  OutputIterator insert(const Vector_3& normal1,
                        const Vector_3& normal2, Vertex_handle vertex2,
                        OutputIterator oi);

  /*! Insert a great arc whose angle is less than Pi and is represented by two
   * normals into the GM. Each normal defines an end point of the great arc.
   * \param normal1 represents the source normal.
   * \param vertex1 the handle of the vertex that is the source of the arc
   * \param normal2 represents the target normal.
   * \param vertex2 the handle of the vertex that is the target of the arc
   * \return the handle for the halfedge directed from the endpoint
   * represented by normal1 toward the endpoint represented by normal2
   */
  template <typename Vector_3, typename OutputIterator>
  OutputIterator insert(const Vector_3& normal1, Vertex_handle vertex1,
                        const Vector_3& normal2, Vertex_handle vertex2,
                        OutputIterator oi);
};

//! \brief inserts a great arc represented by two normals into the GM.
template <typename Vector_3, typename OutputIterator>
OutputIterator Base_gaussian_map_builder::insert(const Vector_3& normal1, const Vector_3& normal2, OutputIterator oi) {
  using Make_x_monotone_result = std::variant<Point, X_monotone_curve>;
  const auto& traits = *(m_gm.geometry_traits());
  std::list<Make_x_monotone_result> x_objects;
  make_x_monotone(normal1, normal2, std::back_inserter(x_objects), traits);
  const auto* xc1 = std::get_if<X_monotone_curve>(&(x_objects.front()));
  Halfedge_handle he = m_gm.insert_in_face_interior(*xc1, m_gm.faces_begin());
  if (! xc1->is_directed_right()) he = he->twin();
  *oi++ = he;
  if (x_objects.size() == 1) return oi;

  const auto* xc2 = std::get_if<X_monotone_curve>(&(x_objects.back()));
  *oi++ = (xc2->is_directed_right()) ?
    m_gm.insert_from_left_vertex(*xc2, he->target()) : m_gm.insert_from_right_vertex(*xc2, he->target());
  return oi;
}

//! \brief inserts a great arc represented by two normals starting at a vertex.
template <typename Vector_3, typename OutputIterator>
OutputIterator Base_gaussian_map_builder::insert(const Vector_3& normal1, Vertex_handle vertex1,
                                                 const Vector_3& normal2, OutputIterator oi) {
  using Make_x_monotone_result = std::variant<Point, X_monotone_curve>;
  const auto& traits = *(m_gm.geometry_traits());
  std::list<Make_x_monotone_result> x_objects;
  make_x_monotone(normal1, normal2, std::back_inserter(x_objects), traits);
  const auto* xc1 = std::get_if<X_monotone_curve>(&(x_objects.front()));
  Halfedge_handle he = (xc1->is_directed_right()) ?
    m_gm.insert_from_left_vertex(*xc1, vertex1) : m_gm.insert_from_right_vertex(*xc1, vertex1);
  *oi++ = he;
  if (x_objects.size() == 1) return oi;

  const auto* xc2 = std::get_if<X_monotone_curve>(&(x_objects.back()));
  *oi++ = (xc2->is_directed_right()) ?
    m_gm.insert_from_left_vertex(*xc2, he->target()) : m_gm.insert_from_right_vertex(*xc2, he->target());
  return oi;
}

//! \brief insert a great arc represented by two normals ending at a vertex.
template <typename Vector_3, typename OutputIterator>
OutputIterator Base_gaussian_map_builder::insert(const Vector_3& normal1, const Vector_3& normal2,
                                                 Vertex_handle vertex2, OutputIterator oi) {
  using Make_x_monotone_result = std::variant<Point, X_monotone_curve>;
  const auto& traits = *(m_gm.geometry_traits());
  std::list<Make_x_monotone_result> x_objects;
  make_x_monotone(normal1, normal2, std::back_inserter(x_objects), traits);
  const auto* xc1 = std::get_if<X_monotone_curve>(&(x_objects.front()));
  if (x_objects.size() == 1) {
    Halfedge_handle he = (xc1->is_directed_right()) ?
      m_gm.insert_from_right_vertex(*xc1, vertex2) : m_gm.insert_from_left_vertex(*xc1, vertex2);
    *oi++ = he->twin();
    return oi;
  }

  const auto* xc2 = std::get_if<X_monotone_curve>(&(x_objects.back()));
  Halfedge_handle he2 = (xc2->is_directed_right()) ?
    m_gm.insert_from_right_vertex(*xc2, vertex2) : m_gm.insert_from_left_vertex(*xc2, vertex2);
  he2 = he2->twin();
  Halfedge_handle he1 = (xc1->is_directed_right()) ?
    m_gm.insert_from_right_vertex(*xc1, he2->source()) : m_gm.insert_from_left_vertex(*xc1, he2->source());
  he1 = he1->twin();
  *oi++ = he1;
  *oi++ = he2;
  return oi;
}

//! \brief inserts a great arc represented by two normals between vertices.
template <typename Vector_3,  typename OutputIterator>
OutputIterator Base_gaussian_map_builder::insert(const Vector_3& normal1, Vertex_handle vertex1,
                                                 const Vector_3& normal2, Vertex_handle vertex2, OutputIterator oi) {
  using Make_x_monotone_result = std::variant<Point, X_monotone_curve>;
  const auto& traits = *(m_gm.geometry_traits());
  std::list<Make_x_monotone_result> x_objects;
  make_x_monotone(normal1, normal2, std::back_inserter(x_objects), traits);
  const auto* xc1 = std::get_if<X_monotone_curve>(&(x_objects.front()));
  if (x_objects.size() == 1) {
    *oi++ = m_gm.insert_at_vertices(*xc1, vertex1, vertex2);
    return oi;
  }

  Halfedge_handle he = (xc1->is_directed_right()) ?
    m_gm.insert_from_left_vertex(*xc1, vertex1) : m_gm.insert_from_right_vertex(*xc1, vertex1);
  *oi++ = he;
  const auto* xc2 = std::get_if<X_monotone_curve>(&(x_objects.back()));
  *oi++ = m_gm.insert_at_vertices(*xc2, he->target(), vertex2);
  return oi;
}

#endif
