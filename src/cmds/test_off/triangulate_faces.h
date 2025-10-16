#ifndef TRIANGULATE_FACES_H
#define TRIANGULATE_FACES_H

#include <list>
#include <unordered_map>

#include <boost/range/value_type.hpp>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Triangulation_data_structure_2.h>

// static std::size_t indent = 0;

/*! Triangulate a face of size at least 5.
 */
template <typename PolygonMesh, typename NamedParameters>
void triangulate_quad(typename boost::graph_traits<PolygonMesh>::face_descriptor fd, PolygonMesh& mesh,
                      const NamedParameters& np) {
  // std::cout << std::string(indent, ' ') << "triangulate_quad()\n";
  using VPM = typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type;
  VPM vpm = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, CGAL::internal_np::vertex_point),
                                               get_property_map(CGAL::vertex_point, mesh));
  auto h0 = CGAL::halfedge(fd, mesh);
  auto h1 = CGAL::next(h0, mesh);
  auto h2 = CGAL::next(h1, mesh);
  auto h3 = CGAL::next(h2, mesh);
  const auto& p0 = get(vpm, CGAL::target(h0, mesh));
  const auto& p1 = get(vpm, CGAL::target(h1, mesh));
  const auto& p2 = get(vpm, CGAL::target(h2, mesh));
  const auto& p3 = get(vpm, CGAL::target(h3, mesh));

  /* Chooses the diagonal that will split the quad in two triangles that
   * maximize the scalar product of of the unnormalized normals of the two
   * triangles.
   * The lengths of the unnormalized normals (computed using
   * cross-products of two vectors) are proportional to the area of the
   * triangles. Maximize the scalar product of the two normals will avoid
   * skinny triangles, and will also take into account the cosine of the
   * angle between the two normals.
   * In particular, if the two triangles are oriented in different
   * directions, the scalar product will be negative.
   */
  auto p21 = p2 - p1;
  auto p12 = p1 - p2;
  auto p32 = p3 - p2;
  auto p03 = p0 - p3;
  auto p30 = p3 - p0;
  auto p10 = p1 - p0;
  auto p1p3 = CGAL::cross_product(p21, p32) * CGAL::cross_product(p03, p10);
  auto p0p2 = CGAL::cross_product(p10, p12) * CGAL::cross_product(p32, p30);
  if (p0p2 > p1p3) CGAL::Euler::split_face(h0, h2, mesh);
  else CGAL::Euler::split_face(h1, h3, mesh);
}

template <typename PolygonMesh>
struct Face_info {
  typename boost::graph_traits<PolygonMesh>::halfedge_descriptor e[3];
  int m_nesting_level;
  bool in_domain() { return m_nesting_level % 2 == 1; }
  // bool is_external;
};

/*! Construct a constrained triangulation
 */
template <typename PolygonMesh, typename Triangulation, typename NamedParameters>
bool construct_triangulation(typename boost::graph_traits<PolygonMesh>::face_descriptor fd,
                             PolygonMesh& mesh, Triangulation& tri, const NamedParameters& np) {
  // std::cout << std::string(indent, ' ') << "construct_triangulation()\n";
  using Graph_traits = boost::graph_traits<PolygonMesh>;
  using vertex_descriptor = typename Graph_traits::vertex_descriptor;
  using halfedge_descriptor = typename Graph_traits::halfedge_descriptor;
  using face_descriptor = typename Graph_traits::face_descriptor;
  using VPM = typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type;
  using Point = typename boost::property_traits<VPM>::value_type;
  VPM vpm = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, CGAL::internal_np::vertex_point),
                                               get_property_map(CGAL::vertex_point, mesh));

  using Vertex_handle = typename Triangulation::Vertex_handle;
  std::unordered_map<Vertex_handle, bool> processed;

  bool res(true);
  auto hd = halfedge(fd, mesh);
  auto range = CGAL::halfedges_around_face(hd, mesh);
  auto range_begin = range.begin();
  auto first_hd = *range_begin;
  auto first_vd = CGAL::target(first_hd, mesh);

  const Point& p = get(vpm, first_vd);
  auto start = tri.insert(p);
  processed[start] = true;

  // std::cout << std::string(indent, ' ') << p << std::endl;

  start->info() = first_hd;
  auto prev = start;

  halfedge_descriptor null_hd = boost::graph_traits<PolygonMesh>::null_halfedge();
  using Range = decltype(range);
  ++range_begin;
  auto skip_first_range = Range(range_begin, range.end());
  for (auto hd : skip_first_range) {
    auto vd = CGAL::target(hd, mesh);

    const Point& p = get(vpm, vd);
    auto next = tri.insert(p);
    auto [it, inserted] = processed.insert({ next, true });

    // std::cout << std::string(indent, ' ') << p << std::endl;

    // next->info() = (next->info() == null_hd) ? hd : null_hd;
    // if (next->info() == null_hd) { ... }
    next->info() = hd;
    if (! inserted) {
      res = false;
      next->info() = null_hd;
      // std::cout << std::string(indent, ' ') << "XXXXX set res to false\n";
    }

    tri.insert_constraint(prev, next);
    prev = next;
  }
  tri.insert_constraint(prev, start);

  processed.clear();
  return res;
}

/*! Mark facets in a triangulation that are inside the domain bounded by
 * the polygon.
 * \param tri (in/out) the triangulation.
 */
template <typename Triangulation>
void mark_domains(Triangulation& tri, typename Triangulation::Face_handle start, int index,
                  std::list<typename Triangulation::Edge>& border) {
  if (start->info().m_nesting_level != -1) return;
  std::list<typename Triangulation::Face_handle> queue;
  queue.push_back(start);
  while (! queue.empty()) {
    auto fh = queue.front();
    queue.pop_front();
    if (fh->info().m_nesting_level == -1) {
      fh->info().m_nesting_level = index;
      for (auto i = 0; i < 3; ++i) {
        typename Triangulation::Edge e(fh,i);
        auto n = fh->neighbor(i);
        if (n->info().m_nesting_level == -1) {
          if (tri.is_constrained(e)) border.push_back(e);
          else queue.push_back(n);
        }
      }
    }
  }
}

  /*! \brief marks facets in a triangulation that are inside the domain.
   * Explores set of facets connected with non constrained edges,
   * and attribute to each such set a nesting level.
   * We start from facets incident to the infinite vertex, with a nesting
   * level of 0. Then we recursively consider the non-explored facets incident
   * to constrained edges bounding the former set and increase the nesting
   * level by 1.
   * Facets in the domain are those with an odd nesting level.
   */
template <typename Triangulation>
void mark_domains(Triangulation& tri) {
  for (auto it = tri.all_faces_begin(); it != tri.all_faces_end(); ++it)
    it->info().m_nesting_level = -1;

  std::list<typename Triangulation::Edge> border;
  mark_domains(tri, tri.infinite_face(), 0, border);
  while (! border.empty()) {
    auto e = border.front();
    border.pop_front();
    auto n = e.first->neighbor(e.second);
    if (n->info().m_nesting_level == -1)
      mark_domains(tri, n, e.first->info().m_nesting_level+1, border);
  }
}

/*! Mark the halfedges on the boundary of a facet as border halfedges.
 */
template <typename PolygonMesh>
void make_hole(typename boost::graph_traits<PolygonMesh>::halfedge_descriptor first_hd, PolygonMesh& mesh) {
  //we are not using Euler::make_hole because it has a precondition
  //that the hole is not made on the boundary of the mesh
  //here we allow making a hole on the boundary, and the pair(s) of
  //halfedges that become border-border are fixed by the connectivity
  //setting made in operator()
  CGAL_assertion(! CGAL::is_border(first_hd, mesh));
  auto fd = CGAL::face(first_hd, mesh);
  for (auto hd : CGAL::halfedges_around_face(first_hd, mesh))
    CGAL::internal::set_border(hd, mesh);
  remove_face(fd, mesh);
}

//!
template <typename FaceHandle>
bool is_external(FaceHandle fh) { return ! fh->info().in_domain(); }

/*! Triangulate a face of size at least 5.
 */
template <typename PolygonMesh, typename Vector_3, typename NamedParameters>
bool triangulate_face(typename boost::graph_traits<PolygonMesh>::face_descriptor fd, PolygonMesh& mesh,
                      const Vector_3 normal, const NamedParameters& np) {
  // std::cout << std::string(indent, ' ') << "triangulate_face(" << normal << ")\n";
  using Graph_traits = boost::graph_traits<PolygonMesh>;
  using vertex_descriptor = typename Graph_traits::vertex_descriptor;
  using halfedge_descriptor = typename Graph_traits::halfedge_descriptor;
  using face_descriptor = typename Graph_traits::face_descriptor;
  using VPM = typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type;
  using Point = typename boost::property_traits<VPM>::value_type;
  using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
  using Traits = CGAL::Projection_traits_3<Kernel>;
  using Vi = typename Graph_traits::halfedge_descriptor;
  using Vb = CGAL::Triangulation_vertex_base_with_info_2<Vi, Traits>;
  using Fi = Face_info<PolygonMesh>;
  using Fbi = CGAL::Triangulation_face_base_with_info_2<Fi, Traits>;
  using Fb = CGAL::Constrained_triangulation_face_base_2<Traits, Fbi>;
  using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
  // typedef CGAL::No_intersection_tag                             Itag;
  using Itag = CGAL::Exact_predicates_tag;
  using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, Tds, Itag>;
  Traits cdt_traits(normal);
  CDT cdt(cdt_traits);
  auto gp = construct_triangulation(fd, mesh, cdt, np);
  mark_domains(cdt);

  // If the facet is not in general position, simply split along one of the
  // edges induced by the triangulation, and recursively triangulate the two
  // resulting facets.
  if (! gp) {
    // std::cout << std::string(indent, ' ') << "triangulate_face() degenerate\n";
    halfedge_descriptor null_hd = boost::graph_traits<PolygonMesh>::null_halfedge();
    halfedge_descriptor hd = null_hd;
    std::size_t i = 0;
    std::size_t j = 0;
    for (auto eit = cdt.finite_edges_begin(); eit != cdt.finite_edges_end(); ++eit) {
      ++i;
      if (cdt.is_constrained(*eit)) continue;

      auto fh = eit->first;
      auto index = eit->second;
      auto fh_opp = fh->neighbor(index);
      if (is_external(fh) || is_external(fh_opp)) continue;
      ++j;

      const auto vh_a = fh->vertex(cdt.cw(index));
      const auto vh_b = fh->vertex(cdt.ccw(index));
      if ((vh_a->info() == null_hd) || (vh_b->info() == null_hd)) continue;

      hd = CGAL::Euler::split_face(vh_a->info(), vh_b->info(), mesh);
      // std::cout << std::string(indent, ' ') << "XXXX 1 triangulate_face(" << normal << ")\n";
      break;
    }
    if (hd == boost::graph_traits<PolygonMesh>::null_halfedge()) {
      std::cerr << "Couldn't find a diagonal!\n";
      return false;
    }
    cdt.clear();

    auto hd_opp = CGAL::opposite(hd, mesh);
    auto fd1 = CGAL::face(hd, mesh);
    auto fd2 = CGAL::face(hd_opp, mesh);

    auto size = CGAL::halfedges_around_face(hd, mesh).size();
    auto size_opp = CGAL::halfedges_around_face(hd_opp, mesh).size();

    bool rc1, rc2;
    if (size == 4) triangulate_quad(fd1, mesh, np);
    else if (size > 4) {
      // indent += 2;
      // std::cout << std::string(indent, ' ') << "XXXX 2 normal: " << normal << std::endl;
      rc1 = triangulate_face(fd1, mesh, normal, np);
      // std::cout << std::string(indent, ' ') << "XXXX 3 normal: " << normal << std::endl;
      // indent -= 2;
    }

    if (size_opp == 4) triangulate_quad(fd2, mesh, np);
    else if (size_opp > 4) {
      // indent += 2;
      rc2 = triangulate_face(fd2, mesh, normal, np);
      // indent -= 2;
    }

    // std::cout << std::string(indent, ' ') << "triangulate_face() degenerate 2 end\n";
    return rc1 || rc2;
  }

  // std::cout << std::string(indent, ' ') << "triangulate_face() non-degenerate\n";
  auto hd = CGAL::halfedge(fd, mesh);
  make_hole(hd, mesh);

  halfedge_descriptor null_hd;
  for (auto eit = cdt.finite_edges_begin(), end = cdt.finite_edges_end(); eit != end; ++eit) {
    auto fh = eit->first;
    auto index = eit->second;
    auto opposite_fh = fh->neighbor(index);
    auto opposite_index = opposite_fh->index(fh);

    const auto vh_a = fh->vertex(cdt.cw(index));
    const auto vh_b = fh->vertex(cdt.ccw(index));
    if ((vh_a->info() == null_hd) || (vh_b->info() == null_hd)) continue;

    //not both fh are external and edge is not constrained
    if ( ! (is_external(fh) && is_external(opposite_fh)) && ! cdt.is_constrained(*eit)) {
      // strictly internal edge
      auto hd_new = CGAL::halfedge(add_edge(mesh), mesh);
      auto hd_opp_new = CGAL::opposite(hd_new, mesh);
      fh->info().e[index] = hd_new;
      opposite_fh->info().e[opposite_index] = hd_opp_new;
      auto hd_a = vh_a->info();
      CGAL::set_target(hd_new, CGAL::target(hd_a, mesh), mesh);
      auto hd_b = vh_b->info();
      CGAL::set_target(hd_opp_new, CGAL::target(hd_b, mesh), mesh);
    }

    if (cdt.is_constrained(*eit)) {
      //edge is constrained
      if (! is_external(fh)) fh->info().e[index] = vh_a->info();
      if (! is_external(opposite_fh))
        opposite_fh->info().e[opposite_index] = vh_b->info();
    }
  }

  for (auto fit = cdt.finite_faces_begin(), end = cdt.finite_faces_end(); fit != end; ++fit) {
    if (is_external(fit)) continue;

    halfedge_descriptor hd0 = fit->info().e[0];
    halfedge_descriptor hd1 = fit->info().e[1];
    halfedge_descriptor hd2 = fit->info().e[2];
    if ((hd0 == null_hd) || (hd1 == null_hd) || (hd2 == null_hd)) {
      std::cerr << "Ignoring a self-intersecting facet\n";
      continue;
    }
    CGAL_assertion(hd0 != null_hd);
    CGAL_assertion(hd1 != null_hd);
    CGAL_assertion(hd2 != null_hd);

    CGAL::set_next(hd0, hd1, mesh);
    CGAL::set_next(hd1, hd2, mesh);
    CGAL::set_next(hd2, hd0, mesh);

    // std::cout << std::string(indent, ' ') << "XXXX before triangulate_face(" << normal << ")\n";
    CGAL::Euler::fill_hole(hd0, mesh);
    // std::cout << std::string(indent, ' ') << "XXXX after triangulate_face(" << normal << ")\n";
  }
  cdt.clear();
  // std::cout << std::string(indent, ' ') << "triangulate_face() non-degenerate end\n";
  return true;
}

/*! Triangulate a mesh
 */
template <typename PolygonMesh, typename Map, typename NamedParameters = CGAL::parameters::Default_named_parameters>
bool triangulate_faces(PolygonMesh& mesh, const Map& normals,
                       const NamedParameters& np = CGAL::parameters::default_values()) {
  // There is a bug in CGAL related to the triangulation of the facets of
  // a face graph. Thus, we use a workaround.
  using Graph_traits = boost::graph_traits<PolygonMesh>;
  using face_descriptor = typename Graph_traits::face_descriptor;

  // Store facet handles, because the collection of the mesh facets
  // is modified during the loop, which invalidates the facet range.
  auto face_range = CGAL::faces(mesh);
  std::vector<face_descriptor> facets;
  facets.reserve(std::distance(boost::begin(face_range), boost::end(face_range)));

  // Filter out triangular faces
  for (auto fd : face_range) {
    auto hd = halfedge(fd, mesh);
    auto size = CGAL::halfedges_around_face(hd, mesh).size();
    if (size > 3) facets.push_back(fd);
  }

  // Iterates on the vector of face descriptors
  auto res = true;
  for (auto fd : facets) {
    auto hd = halfedge(fd, mesh);
    auto size = CGAL::halfedges_around_face(hd, mesh).size();
    if (size == 4) {
      triangulate_quad(fd, mesh, np);
      continue;
    }
    if (! triangulate_face(fd, mesh, get(normals, fd), np)) res = false;
  }
  return res;
}

#endif
