#ifndef MATCH_FACES_H
#define MATCH_FACES_H

#include <vector>
#include <utility>
#include <algorithm>

#include <boost/container/small_vector.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/dynamic_bitset.hpp>

#include <CGAL/assertions.h>
#include <CGAL/boost/graph/iterator.h>
#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/properties.h>
#include <CGAL/boost/graph/named_params_helper.h>
#include <CGAL/Named_function_parameters.h>

// #include <CGAL/Polygon_mesh_processing/border.h>
// #include <CGAL/utils_classes.h>

namespace CGAL {
namespace Polygon_mesh_processing {

////////

template< typename PolygonMesh1,
          typename PolygonMesh2,
          typename FacePairOutputIterator,
          typename FaceOutputIterator1,
          typename FaceOutputIterator2,
          typename NamedParameters1 = parameters::Default_named_parameters,
          typename NamedParameters2 = parameters::Default_named_parameters >
void my_match_faces(const PolygonMesh1& m1,
                 const PolygonMesh2& m2,
                 FacePairOutputIterator common,
                 FaceOutputIterator1 m1_only,
                 FaceOutputIterator2 m2_only,
                 const NamedParameters1& np1 = parameters::default_values(),
                 const NamedParameters2& np2 = parameters::default_values())
{
  typedef typename GetVertexPointMap<PolygonMesh1, NamedParameters1>::const_type            VPMap1;
  typedef typename GetVertexPointMap<PolygonMesh2, NamedParameters2>::const_type            VPMap2;
  typedef typename GetInitializedVertexIndexMap<PolygonMesh1, NamedParameters1>::const_type VIMap1;
  typedef typename GetInitializedVertexIndexMap<PolygonMesh2, NamedParameters2>::const_type VIMap2;
  typedef typename boost::property_traits<VPMap2>::value_type                               Point_3;
  typedef typename boost::graph_traits<PolygonMesh1>::face_descriptor                       face_descriptor_1;

  using parameters::choose_parameter;
  using parameters::get_parameter;

  const VPMap1 vpm1 = choose_parameter(get_parameter(np1, internal_np::vertex_point),
                                       get_const_property_map(vertex_point, m1));
  const VPMap2 vpm2 = choose_parameter(get_parameter(np2, internal_np::vertex_point),
                                       get_const_property_map(vertex_point, m2));
  static_assert(std::is_same<typename boost::property_traits<VPMap1>::value_type,
                             typename boost::property_traits<VPMap2>::value_type>::value,
                            "Both vertex point maps must have the same point type.");

  const VIMap1 vim1 = get_initialized_vertex_index_map(m1, np1);
  const VIMap2 vim2 = get_initialized_vertex_index_map(m2, np2);

  std::map<Point_3, std::size_t> point_id_map;

  std::vector<std::size_t> m1_vertex_id(num_vertices(m1), -1);
  std::vector<std::size_t> m2_vertex_id(num_vertices(m2), -1);
  boost::dynamic_bitset<> shared_vertices(m1_vertex_id.size() + m2_vertex_id.size());

  //iterate both meshes to set ids of all points, and set vertex/point_id maps.
  std::size_t id = 0;
  for(auto v : vertices(m1))
  {
    const typename boost::property_traits<VPMap1>::reference p = get(vpm1, v);
    auto res = point_id_map.emplace(p, id);
    if(res.second)
      ++id;
    m1_vertex_id[get(vim1, v)] = res.first->second;
  }
  for(auto v : vertices(m2))
  {
    const typename boost::property_traits<VPMap2>::reference p = get(vpm2, v);
    auto res = point_id_map.emplace(p, id);
    if(res.second)
      ++id;
    else
      shared_vertices.set(res.first->second);
    m2_vertex_id[get(vim2, v)] = res.first->second;
  }

  //fill a set with the "faces point-ids" of m1 and then iterate faces of m2 to compare.
  std::map<boost::container::small_vector<std::size_t, 4>, face_descriptor_1> m1_faces_map;
  for(auto f : faces(m1))
  {
    bool all_shared = true;
    boost::container::small_vector<std::size_t, 4> ids;
    for(auto hd : CGAL::halfedges_around_face(halfedge(f, m1), m1))
    {
      if (CGAL::face(hd, m1) == CGAL::face(CGAL::opposite(hd, m1), m1)) {
        continue;
      }
      auto v = CGAL::target(hd, m1);
      std::size_t vid = m1_vertex_id[get(vim1, v)];
      ids.push_back(vid);
      if(!shared_vertices.test(vid))
      {
        all_shared = false;
        break;
      }
    }
    if(all_shared)
    {
      internal::rearrange_face_ids(ids);
      m1_faces_map.emplace(ids, f);
    }
    else {
      *m1_only++ = f;
    }
  }
  for(auto f : faces(m2))
  {
    boost::container::small_vector<std::size_t, 4> ids;
    bool all_shared = true;
    for(auto hd : CGAL::halfedges_around_face(halfedge(f, m2), m2))
    {
      if (CGAL::face(hd, m2) == CGAL::face(CGAL::opposite(hd, m2), m2)) {
        continue;
      }
      auto v = CGAL::target(hd, m2);
      std::size_t vid = m2_vertex_id[get(vim2, v)];
      ids.push_back(vid);
      if(!shared_vertices.test(vid))
      {
        all_shared = false;
        break;
      }
    }
    if(all_shared)
    {
      internal::rearrange_face_ids(ids);
      auto it = m1_faces_map.find(ids);
      if(it != m1_faces_map.end())
      {
        *common++ = std::make_pair(it->second, f);
        m1_faces_map.erase(it);
      }
      else
      {
        // std::cout << "XXXXXXXX 3\n";
        // for (auto id : ids) std::cout << id << " ";
        // std::cout << std::endl;

        *m2_only++ = f;
      }
    }
    else {
      *m2_only++ = f;
    }
  }
  //all shared faces have been removed from the map, so all that remains must go in m1_only
  for(const auto& it : m1_faces_map)
  {
    // std::cout << "XXXXXXXX 4\n";
    // const auto& ids = it.first;
    // for (auto id : ids) std::cout << id << " ";
    // std::cout << std::endl;

    *m1_only++ = it.second;
  }
}

} // namespace Polygon_mesh_processing
} // namespace CGAL

#endif
