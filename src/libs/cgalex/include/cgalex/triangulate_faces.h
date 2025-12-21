#ifndef TRIANGULATE_FACES_H
#define TRIANGULATE_FACES_H

#include <unordered_set>
#include <unordered_map>
#include <stack>

#include <boost/graph/graph_traits.hpp>
#include <boost/range/value_type.hpp>

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Projection_traits_3.h>
#include <CGAL/Named_function_parameters.h>
#include <CGAL/Triangulation_data_structure_2.h>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

template <typename PolygonMesh, typename VPM>
class Face_triangulator {
private:
  using Polygon_mesh = PolygonMesh;
  using Graph_traits = boost::graph_traits<Polygon_mesh>;
  using vertex_descriptor = typename Graph_traits::vertex_descriptor;
  using halfedge_descriptor = typename Graph_traits::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<Polygon_mesh>::face_descriptor;

  Polygon_mesh& m_mesh;
  const VPM& m_vpm;

  /*! Construct a constrained triangulation
   */
  template <typename Triangulation>
  void construct_triangulation(face_descriptor fd, Triangulation& tri) {
    // std::cout << "construct_triangulation(" << fd << ")\n";

    using Point = typename boost::property_traits<VPM>::value_type;

    // We need to identify vertices of the triangulation that are mapped to the same point in 2D.
    // Observe that a triangulation vertex is associated with at least one mesh vertex.
    // In a degenerate case, a triangulation vertex can be associated with several mesh vertices.
    // The insertion of a point into the triangulation returns a triangulation vertex handle.
    // Unfortunately, there is no indication whether the point has already been inserted.

    using Vertex_handle = typename Triangulation::Vertex_handle;
    auto null_vh = Vertex_handle();

    // Obtain a handle to the descriptor of the first halfedge that is not degenerate
    auto rep_hd = halfedge(fd, m_mesh);
    auto range = CGAL::halfedges_around_face(rep_hd, m_mesh);
    auto it = range.begin();
    for (; it != range.end(); ++it) {
      auto hd = *it;
      auto ohd = CGAL::opposite(hd, m_mesh);
      auto ofd = CGAL::face(ohd, m_mesh);
      if (fd != ofd) break;
    }

    // Insert the first point into the triangulation
    auto hd = *it;
    auto prev_vd = CGAL::source(hd, m_mesh);
    const Point& prev_p = get(m_vpm, prev_vd);
    auto prev_vh = tri.insert(prev_p);

    while (it != range.end()) {
      // Insert the second point into the triangulation
      auto vd = CGAL::target(hd, m_mesh);
      const Point& p = get(m_vpm, vd);
      auto next_vh = tri.insert(p);
      next_vh->info().insert(hd);

      // Insert the constraints
      // std::cout << prev_vh->point() << ", " << next_vh->point() << std::endl;
      tri.insert_constraint(prev_vh, next_vh);
      prev_vh = next_vh;

      // Advance
      bool degenerate = false;
      for (++it; it != range.end(); ++it) {
        hd = *it;
        auto ohd = CGAL::opposite(hd, m_mesh);
        auto ofd = CGAL::face(ohd, m_mesh);
        if (fd != ofd) break;
        degenerate = true;
      }
      if (! degenerate) continue;
      if (it == range.end()) break;

      auto prev_vd = CGAL::source(hd, m_mesh);
      const Point& prev_p = get(m_vpm, prev_vd);
      prev_vh = tri.insert(prev_p);
    }
    // std::cout << "End construct_triangulation(" << fd << ")\n";
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
  void make_hole(halfedge_descriptor first_hd) {
    //we are not using Euler::make_hole because it has a precondition
    //that the hole is not made on the boundary of the mesh
    //here we allow making a hole on the boundary, and the pair(s) of
    //halfedges that become border-border are fixed by the connectivity
    //setting made in operator()
    CGAL_assertion(! CGAL::is_border(first_hd, m_mesh));
    auto fd = CGAL::face(first_hd, m_mesh);
    for (auto hd : CGAL::halfedges_around_face(first_hd, m_mesh))
      CGAL::internal::set_border(hd, m_mesh);
    remove_face(fd, m_mesh);
  }

  struct Face_info {
    halfedge_descriptor e[3];
    int m_nesting_level;
    bool in_domain() { return m_nesting_level % 2 == 1; }
    // bool is_external;
  };

  //!
  template <typename FaceHandle>
  bool is_external(FaceHandle fh) { return ! fh->info().in_domain(); }

public:
  Face_triangulator(Polygon_mesh& mesh, const VPM& vpm) :
    m_mesh(mesh),
    m_vpm(vpm)
  {}

  //
  void my_remove_face(halfedge_descriptor h) {
    // std::cout << "my_remove_face()\n";
    CGAL_precondition(CGAL::is_valid_halfedge_descriptor(h, m_mesh));
    CGAL_precondition(! CGAL::is_border(h, m_mesh));

    auto f = CGAL::face(h, m_mesh);

    // Advance around the face until a non-border opposite halfedge is encountered.
    halfedge_descriptor end = h;
    do {
      halfedge_descriptor nh = CGAL::next(h, m_mesh);
      halfedge_descriptor oh = CGAL::opposite(h, m_mesh);
      face_descriptor of = CGAL::face(oh, m_mesh);
      if (! CGAL::is_border(oh, m_mesh) && (f != of)) break;
      h = nh;
    } while (h != end);

    // If a non-border opposite halfedge is not encountered, remove everything.
    if ((h == end) && CGAL::is_border(CGAL::opposite(h, m_mesh), m_mesh)) {
      do {
        halfedge_descriptor nh = CGAL::next(h, m_mesh);
        CGAL::remove_vertex(CGAL::target(h, m_mesh), m_mesh);
        CGAL::remove_edge(CGAL::edge(h, m_mesh), m_mesh);
        h = nh;
      } while (h != end);
      CGAL::remove_face(f, m_mesh);
      return;
    }

    // Make sure that h is not a border halfedge!
    CGAL_precondition(! CGAL::is_border(CGAL::opposite(h, m_mesh), m_mesh) &&
                      (f != CGAL::face(CGAL::opposite(h, m_mesh), m_mesh)));
    end = h;
    do {
      halfedge_descriptor nh = CGAL::next(h, m_mesh);
      halfedge_descriptor oh = CGAL::opposite(h, m_mesh);
      bool h_border = CGAL::is_border(oh, m_mesh);

      CGAL::internal::set_border(h, m_mesh);
      // CGAL::set_face(h, boost::graph_traits<PolygonMesh>::null_face(), m_mesh);

      if (! h_border) {
        h = nh;
        continue;
      }

      // Remove the edge
      // Handle the source vertex
      halfedge_descriptor ph = CGAL::prev(h, m_mesh);
      if (ph != oh) {
        set_halfedge(target(oh, m_mesh), ph, m_mesh);
        set_next(ph, CGAL::next(oh, m_mesh), m_mesh);
      }
      else {
        CGAL::remove_vertex(CGAL::target(oh, m_mesh), m_mesh);
      }

      // Handle the target vertex
      if (nh != oh) {
        halfedge_descriptor onh = CGAL::opposite(nh, m_mesh);
        set_halfedge(target(h, m_mesh), onh, m_mesh);
        halfedge_descriptor poh = CGAL::prev(oh, m_mesh);
        set_next(poh, nh, m_mesh);
        h = nh;
      }
      else {
        CGAL::remove_vertex(CGAL::target(h, m_mesh), m_mesh);
        h = CGAL::next(nh, m_mesh);
      }

      CGAL::remove_edge(CGAL::edge(oh, m_mesh), m_mesh);
    } while (h != end);

    CGAL::remove_face(f, m_mesh);
    // std::cout << "end my_remove_face()\n";
  }

  /*! Triangulate a face of size at least 5.
   */
  void triangulate_quad(face_descriptor fd) {
    auto h0 = CGAL::halfedge(fd, m_mesh);
    auto h1 = CGAL::next(h0, m_mesh);
    auto h2 = CGAL::next(h1, m_mesh);
    auto h3 = CGAL::next(h2, m_mesh);
    const auto& p0 = get(m_vpm, CGAL::target(h0, m_mesh));
    const auto& p1 = get(m_vpm, CGAL::target(h1, m_mesh));
    const auto& p2 = get(m_vpm, CGAL::target(h2, m_mesh));
    const auto& p3 = get(m_vpm, CGAL::target(h3, m_mesh));

    /* Choose the diagonal that will split the quad in two triangles that
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

    if (p0p2 > p1p3) CGAL::Euler::split_face(h0, h2, m_mesh);
    else CGAL::Euler::split_face(h1, h3, m_mesh);
  }

  /*! Triangulate a face of size at least 5.
   */
  template <typename Vector_3>
  bool triangulate_face(face_descriptor fd, const Vector_3 normal) {
    // std::cout << "triangulate_face(" << fd << ", " << normal << ")\n";

    using Point = typename boost::property_traits<VPM>::value_type;
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Traits = CGAL::Projection_traits_3<Kernel>;
    using Vi = std::unordered_set<halfedge_descriptor>;
    using Vb = CGAL::Triangulation_vertex_base_with_info_2<Vi, Traits>;
    using Fi = Face_info;
    using Fbi = CGAL::Triangulation_face_base_with_info_2<Fi, Traits>;
    using Fb = CGAL::Constrained_triangulation_face_base_2<Traits, Fbi>;
    using Tds = CGAL::Triangulation_data_structure_2<Vb, Fb>;
    // using Itag = CGAL::No_intersection_tag;
    using Itag = CGAL::Exact_predicates_tag;
    using CDT = CGAL::Constrained_Delaunay_triangulation_2<Traits, Tds, Itag>;
    Traits cdt_traits(normal);
    CDT cdt(cdt_traits);
    construct_triangulation(fd, cdt);

    mark_domains(cdt);

    auto hd = CGAL::halfedge(fd, m_mesh);
    halfedge_descriptor null_hd;

    // for (auto vit = cdt.vertices_begin(); vit != cdt.vertices_end(); ++vit) {
    //   const auto& hds = vit->info();
    //   std::cout << vit->point() << ", " << hds.size() << std::endl;
    // }

    // std::cout << m_mesh.number_of_vertices() << "," << m_mesh.number_of_halfedges() << ","
    //           << m_mesh.number_of_faces() << std::endl;

    // Remove the face along with the internal edges and vertices.
    // CGAL::Euler::remove_face(hd, m_mesh);
    my_remove_face(hd);

    // std::cout << m_mesh.number_of_vertices() << "," << m_mesh.number_of_halfedges() << ","
    //           << m_mesh.number_of_faces() << std::endl;

    for (auto eit = cdt.finite_edges_begin(), end = cdt.finite_edges_end(); eit != end; ++eit) {
      auto fh = eit->first;
      auto index = eit->second;
      auto opposite_fh = fh->neighbor(index);
      auto opposite_index = opposite_fh->index(fh);

      const auto vh_a = fh->vertex(cdt.cw(index));
      const auto vh_b = fh->vertex(cdt.ccw(index));

      // none of fh and fh_opposite are external and edge is not constrained
      if (! is_external(fh) && ! is_external(opposite_fh) && ! cdt.is_constrained(*eit)) {
        // strictly internal edge
        const auto& hd_as = vh_a->info();
        const auto& hd_bs = vh_b->info();
        // std::cout << "Internal: " << hd_as.size() << ", " << hd_bs.size() << std::endl;

        // const auto& pa = vh_a->point();
        // const auto& pb = vh_b->point();
        // std::cout << pa << " " << pb << std::endl;

        const auto& hdas = vh_a->info();
        auto ita = hdas.begin();
        auto hd_a = *ita;
        // auto next_hd_a = CGAL::next(hd_a, m_mesh);
        // auto p_vd = CGAL::source(hd_a, m_mesh);
        // const Point& p = get(m_vpm, p_vd);
        // auto q_vd = CGAL::target(next_hd_a, m_mesh);
        // const Point& q = get(m_vpm, q_vd);
        // std::cout << "p: " << p << std::endl;
        // std::cout << "q: " << q << std::endl;
        for (++ita; ita != hdas.end(); ++ita) {
          if (true) {
            std::cout << "XXXXXXXX 1\n";
            break;
          }
        }

        const auto& hdbs = vh_b->info();
        auto itb = hdbs.begin();
        auto hd_b = *itb;
        // auto next_hd_b = CGAL::next(hd_b, m_mesh);
        for (++itb; itb != hdbs.end(); ++itb) {
          if (true) {
            std::cout << "XXXXXXXX 2\n";
            break;
          }
        }

        auto hd_new = CGAL::halfedge(add_edge(m_mesh), m_mesh);
        auto hd_opp_new = CGAL::opposite(hd_new, m_mesh);
        fh->info().e[index] = hd_new;
        opposite_fh->info().e[opposite_index] = hd_opp_new;
        CGAL::set_target(hd_new, CGAL::target(hd_a, m_mesh), m_mesh);
        CGAL::set_target(hd_opp_new, CGAL::target(hd_b, m_mesh), m_mesh);
      }

      if (cdt.is_constrained(*eit)) {
        // edge is constrained
        const auto& hd_as = vh_a->info();
        const auto& hd_bs = vh_b->info();
        // std::cout << "Constrained: " << hd_as.size() << ", " << hd_bs.size() << std::endl;
        // std::cout << "Edge: " << vh_a->point() << "," << vh_b->point() << std::endl;
        if (! is_external(fh)) {
          const auto& hds = vh_a->info();
          auto it = hds.begin();
          auto hd_a = *it;
          for (++it; it != hds.end(); ++it) {
            auto src_vd = CGAL::source(*it, m_mesh);
            const Point& src = get(m_vpm, src_vd);
            if (src == vh_b->point()) {
              hd_a = *it;
              break;
            }
          }
          fh->info().e[index] = hd_a;
        }

        //
        if (! is_external(opposite_fh)) {
          const auto& hds = vh_b->info();
          auto it = hds.begin();
          auto hd_b = *it;
          for (++it; it != hds.end(); ++it) {
            auto src_vd = CGAL::source(*it, m_mesh);
            const Point& src = get(m_vpm, src_vd);
            if (src == vh_a->point()) {
              hd_b = *it;
              break;
            }
          }
          opposite_fh->info().e[opposite_index] = hd_b;
        }
      }
    }

    for (auto fit = cdt.finite_faces_begin(), end = cdt.finite_faces_end(); fit != end; ++fit) {
      if (is_external(fit)) continue;

      // if (cnt == 1) break;

      halfedge_descriptor hd0 = fit->info().e[0];
      halfedge_descriptor hd1 = fit->info().e[1];
      halfedge_descriptor hd2 = fit->info().e[2];

      if ((hd0 == null_hd) || (hd1 == null_hd) || (hd2 == null_hd)) {
        exit(1);
        continue;
      }
      CGAL_assertion(hd0 != null_hd);
      CGAL_assertion(hd1 != null_hd);
      CGAL_assertion(hd2 != null_hd);

      CGAL::set_next(hd0, hd1, m_mesh);
      CGAL::set_next(hd1, hd2, m_mesh);
      CGAL::set_next(hd2, hd0, m_mesh);

      CGAL::Euler::fill_hole(hd0, m_mesh);
    }

    cdt.clear();
    // std::cout << "End triangulate_face(" << fd << ")\n";
    return true;
  }
};

/*! Triangulate a mesh
 */
template <typename PolygonMesh, typename Map, typename NamedParameters = CGAL::parameters::Default_named_parameters>
bool triangulate_faces(PolygonMesh& mesh, const Map& normals, const NamedParameters& np = params::default_values()) {
  // std::cout << "triangulate_faces(()\n";
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

  using VPM = typename CGAL::GetVertexPointMap<PolygonMesh, NamedParameters>::type;
  VPM vpm = CGAL::parameters::choose_parameter(CGAL::parameters::get_parameter(np, CGAL::internal_np::vertex_point),
                                               get_property_map(CGAL::vertex_point, mesh));
  Face_triangulator triangulator(mesh, vpm);

  // Iterates on the vector of face descriptors
  auto res = true;
  for (auto fd : facets) {
    auto hd = halfedge(fd, mesh);
    auto size = CGAL::halfedges_around_face(hd, mesh).size();
    if (size == 4) {
      triangulator.triangulate_quad(fd);
      continue;
    }
    if (! triangulator.triangulate_face(fd, get(normals, fd))) res = false;
  }
  // std::cout << "End triangulate_faces(()\n";
  return res;
}

#endif
