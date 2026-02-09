#ifndef CGALEX_SIMPLE_GAUSSIAN_MAP_BUILDER_H
#define CGALEX_SIMPLE_GAUSSIAN_MAP_BUILDER_H

#include <unordered_map>
#include <list>
#include <algorithm>

#include <boost/graph/graph_traits.hpp>

#include <CGAL/Kernel_traits.h>

#include "cgalex/Base_gaussian_map_builder.h"

template <typename Arrangement_, typename Graph_, typename FaceNormalMap, typename VertexPointMap>
class Simple_gaussian_map_builder : public Base_gaussian_map_builder<Arrangement_> {
public:
  using Arrangement = Arrangement_;
  using Graph = Graph_;
  using Face_normal_map = FaceNormalMap;
  using Vertex_point_map = VertexPointMap;

  using Vertex_handle = typename Arrangement::Vertex_handle;
  using Halfedge_handle = typename Arrangement::Halfedge_handle;
  using Face_handle = typename Arrangement::Face_handle;
  using Graph_traits = typename boost::graph_traits<Graph>;
  using vertex_descriptor = typename Graph_traits::vertex_descriptor;
  using halfedge_descriptor = typename Graph_traits::halfedge_descriptor;
  using face_descriptor = typename Graph_traits::face_descriptor;
  using halfedge_around_target_circulator = CGAL::Halfedge_around_target_circulator<Graph>;
  using Point_3 = typename boost::property_traits<VertexPointMap>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Point_3>::Kernel;
  using Vector_3 = typename Kernel::Vector_3;

  using Base = Base_gaussian_map_builder<Arrangement>;

  /*! Construct
   */
  Simple_gaussian_map_builder(Arrangement& gm, const Face_normal_map& normals, const Vertex_point_map& points);

  /*! Build the Gaussian map of a graph.
   * \param src the graph vertex currently processed
   * \param graph the graph
   * \param first_time true if the invocation to this function is the first
   */
  void operator()(const Graph& graph, vertex_descriptor src, bool first_time = true);

private:
  /*! Rocesses a graph halfedge
   */
  void process(const Graph& g, vertex_descriptor src, halfedge_around_target_circulator hec,
               halfedge_around_target_circulator next_hec, bool first_time);

  // Data memebers
  const Face_normal_map& m_normals;
  const Vertex_point_map& m_points;
  std::unordered_map<vertex_descriptor, bool> m_vprocessed;
  std::unordered_map<halfedge_descriptor, bool> m_hprocessed;
  std::unordered_map<face_descriptor, Vertex_handle> m_vertex_handles;
};

//! \brief constructs
template <typename Arr, typename G, typename Fnm, typename Vpm>
Simple_gaussian_map_builder<Arr, G, Fnm, Vpm>::Simple_gaussian_map_builder(Arrangement& gm,
                                                                           const Face_normal_map& normals,
                                                                           const Vertex_point_map& points) :
  Base(gm),
  m_normals(normals),
  m_points(points)
{}

//! \brief processes a graph halfedge
template <typename Arr, typename G, typename Fnm, typename Vpm>
void Simple_gaussian_map_builder<Arr, G, Fnm, Vpm>::process(const Graph& g, vertex_descriptor src,
                                                            halfedge_around_target_circulator hec,
                                                            halfedge_around_target_circulator next_hec,
                                                            bool first_time) {
  Vertex_handle invalid_vertex;

  const auto& normal1 = get(m_normals, CGAL::face(*hec, g));
  const auto& normal2 = get(m_normals, CGAL::face(*next_hec, g));
  auto trg = CGAL::target(CGAL::opposite(*next_hec, g), g);
  // The arc might be non-x-monotone. In this case, it is broken into 2
  // x-monotone curves. The halfedges of both are obtained.
  std::list<Halfedge_handle> hes;
  if (first_time) {
    this->insert(normal1, normal2, std::back_inserter(hes));
    m_vertex_handles[CGAL::face(*hec, g)] = hes.front()->source();
    m_vertex_handles[CGAL::face(*next_hec, g)] = hes.back()->target();
    return;
  }

  Vertex_handle v1 = m_vertex_handles[CGAL::face(*hec, g)];
  Vertex_handle v2 = m_vertex_handles[CGAL::face(*next_hec, g)];
  if ((v1 != invalid_vertex) && (v2 != invalid_vertex)) {
    this->insert(normal1, v1, normal2, v2, std::back_inserter(hes));
    Face_handle src_face = hes.front()->twin()->face();
    Face_handle trg_face = hes.front()->face();
    src_face->set_data(get(m_points, src));
    trg_face->set_data(get(m_points, trg));
    return;
  }
  if (v1 != invalid_vertex) {
    this->insert(normal1, v1, normal2, std::back_inserter(hes));
    m_vertex_handles[CGAL::face(*next_hec, g)] = hes.back()->target();
    return;
  }
  if (v2 != invalid_vertex) {
    this->insert(normal1, normal2, v2, std::back_inserter(hes));
    m_vertex_handles[CGAL::face(*hec, g)] = hes.front()->source();
    return;
  }
  CGAL_error_msg("Error: in process(); both vertices are invalid!");
}

//! \brief processes a graph vertex recursively constructing the GM
template <typename Arr, typename G, typename Fnm, typename Vpm>
void Simple_gaussian_map_builder<Arr, G, Fnm, Vpm>::operator()(const Graph& g, vertex_descriptor src, bool first_time) {
  CGAL::Halfedge_around_target_circulator<Graph> hec(src, g);
  CGAL_assertion(CGAL::circulator_size(hec) >= 3);

  // If this is not the first invocation, advance the halfedge iterator until
  // its source vertex is processed. It is guaranteed to reach such a halfedge
  // on consecutive invocations.
  if (! first_time) while (! m_vprocessed[CGAL::source(*hec, g)]) ++hec;

  // Traverse the incident halfedges:
  auto next_hec = hec;
  auto begin_hec = next_hec++;
  do {
    if (! m_hprocessed[*next_hec]) {
      process(g, src, hec, next_hec, first_time);
      first_time = false;
      m_hprocessed[*next_hec] = true;
      m_hprocessed[CGAL::opposite(*next_hec, g)] = true;
    }
    hec = next_hec;
    ++next_hec;
  } while (hec != begin_hec);
  m_vprocessed[src] = true;

  // Traverse recursively:
  hec = CGAL::Halfedge_around_target_circulator<Graph>(src, g);
  begin_hec = hec;
  do {
    auto vd = CGAL::target(CGAL::opposite(*hec, g), g);
    if (! m_vprocessed[vd]) operator()(g, vd, false);
  } while (++hec != begin_hec);
}

#endif
