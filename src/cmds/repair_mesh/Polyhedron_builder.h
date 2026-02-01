#ifndef POLYHEDRON_BUILDER_H
#define POLYHEDRON_BUILDER_H

#include <unordered_map>

#include <CGAL/Polyhedron_incremental_builder_3.h>

template <typename Hds, typename GaussianMap>
class Polyhedron_builder : public CGAL::Modifier_base<Hds> {
private:
  using Gaussian_map = GaussianMap;
  const Gaussian_map& m_gm;

public:
  Polyhedron_builder(const Gaussian_map& gm) : m_gm(gm) {}

  void operator()(Hds& hds) {
    using Face_const_handle = typename Gaussian_map::Face_const_handle;
    CGAL::Polyhedron_incremental_builder_3<Hds> B(hds, true);
    std::size_t num_vertices = m_gm.number_of_faces();
    std::size_t num_facets = 0;
    std::unordered_map<Face_const_handle, std::size_t> indices;
    for (auto v : m_gm.vertex_handles()) if (v->degree() > 2) ++num_facets;
    B.begin_surface(num_vertices, num_facets);
    std::size_t index = 0;
    for (auto f : m_gm.face_handles()) {
      B.add_vertex(f->data());
      indices[f] = index++;
    }
    for (auto v : m_gm.vertex_handles()) {
      if (v->degree() == 2) continue;
      B.begin_facet();
      auto first = v->incident_halfedges();
      auto curr = first;
      do B.add_vertex_to_facet(indices[curr->face()]);
      while (--curr != first);
      B.end_facet();
    }
    B.end_surface();
  }
};

#endif
