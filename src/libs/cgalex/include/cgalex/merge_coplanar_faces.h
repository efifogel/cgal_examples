#ifndef MERGE_COPLANAR_FACES_H
#define MERGE_COPLANAR_FACES_H

#include <list>
#include <unordered_map>

#include <boost/graph/graph_traits.hpp>

#include <CGAL/Union_find.h>

/*!
 */
template <typename FaceDescriptor, typename HandleMap>
typename CGAL::Union_find<FaceDescriptor>::handle
uf_get_face_handle(FaceDescriptor f, CGAL::Union_find<FaceDescriptor>& uf_faces, HandleMap& handles) {
  using Face_descriptor = FaceDescriptor;
  using Uf_faces = CGAL::Union_find<Face_descriptor>;
  using Uf_handle = typename Uf_faces::handle;

  auto insert_res = handles.emplace(f, Uf_handle());
  if (insert_res.second) insert_res.first->second = uf_faces.make_set(f);
  return insert_res.first->second;
}

/*!
 */
template <typename FaceDescriptor, typename HandleMap>
bool uf_join_faces(FaceDescriptor f1, FaceDescriptor f2,
                   CGAL::Union_find<FaceDescriptor>& uf_faces,
                   HandleMap& handles) {
  using Face_descriptor = FaceDescriptor;
  using Uf_faces = CGAL::Union_find<Face_descriptor>;
  using Uf_handle = typename Uf_faces::handle;

  Uf_handle h1 = uf_get_face_handle(f1, uf_faces, handles);
  Uf_handle h2 = uf_get_face_handle(f2, uf_faces, handles);
  bool same_set = uf_faces.same_set(h1, h2);
  if (! same_set) uf_faces.unify_sets(h1, h2);
  return ! same_set;
}

/*! Remove internal vertices. In particular, traverse all vertices and remove internal vertices of degree 2.
 * We assume that the faces are not necessarily convex.
 */
template <typename Graph, typename NamedParameters = CGAL::parameters::Default_named_parameters>
void remove_internal_vertices(Graph& graph, const NamedParameters& np = CGAL::parameters::default_values()) {
  for (auto vd : CGAL::vertices(graph)) {
    if (CGAL::degree(vd, graph) != 2) continue;
    auto hd = CGAL::halfedge(vd, graph);
    auto fd = CGAL::face(hd, graph);
    auto ohd = CGAL::opposite(hd, graph);
    auto ofd = CGAL::face(ohd, graph);
    if (fd == ofd) continue;
    CGAL::Euler::join_vertex(CGAL::opposite(hd, graph), graph);
  }
}

/*! Merge coplanar facets.
 */
template <typename Graph, typename Map, typename NamedParameters>
void merge_coplanar_faces(Graph& graph, const Map& normals,
                          const NamedParameters& np = CGAL::parameters::default_values()) {
  namespace parms = CGAL::parameters;

  using Named_parameters = NamedParameters;
  using Graph_traits = boost::graph_traits<Graph>;
  using halfedge_descriptor = typename Graph_traits::halfedge_descriptor;
  using face_descriptor = typename Graph_traits::face_descriptor;
  using Uf_faces = CGAL::Union_find<face_descriptor>;
  using Kernel = typename CGAL::GetGeomTraits<Graph, Named_parameters>::type;
  const auto& param = parms::get_parameter(np, CGAL::internal_np::geom_traits);
  auto kernel = parms::choose_parameter<Kernel>(param);
  Uf_faces uf_faces;
  std::unordered_map<face_descriptor, typename Uf_faces::handle> uf_handles;
  std::list<halfedge_descriptor> hds;
  auto eq = kernel.equal_3_object();

  // Group coplanar facets in subsets
  for (auto ed : CGAL::edges(graph)) {
    auto hd = CGAL::halfedge(ed, graph);
    auto fd = CGAL::face(hd, graph);
    auto oed = CGAL::opposite(hd, graph);
    auto ofd = CGAL::face(oed, graph);
    const auto& dir = get(normals, fd).direction();
    const auto& odir = get(normals, ofd).direction();
    if ((fd != ofd) && eq(dir, odir)) {
      if (uf_join_faces(fd, ofd, uf_faces, uf_handles)) hds.emplace_back(hd);
    }
  }

  for (auto hd : hds) CGAL::Euler::join_face(hd, graph);

  // Traverse all vertices and remove antenas.
  bool done;
  do {
    done = true;
    for (auto vd : CGAL::vertices(graph)) {
      if (CGAL::degree(vd, graph) != 1) continue;
      auto hd = CGAL::halfedge(vd, graph);
      CGAL::Euler::remove_center_vertex(hd, graph);
      done = false;
    }
  } while (! done);

  // Remove internal vertices
  remove_internal_vertices(graph, np);
}

#endif
