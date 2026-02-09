#ifndef CGALEX_BUILD_GAUSSIAN_MAP_H
#define CGALEX_BUILD_GAUSSIAN_MAP_H

#include <unordered_map>

#include <boost/property_map/property_map.hpp>

#include <CGAL/Named_function_parameters.h>

#include "cgalex/Gaussian_map_builder.h"

namespace cp = CGAL::parameters;
namespace ci = CGAL::internal_np;

template <typename Graph, typename FaceNormalMap, typename Arrangement,
          typename NamedParameters = cp::Default_named_parameters>
void build_gaussian_map(Graph& g, const FaceNormalMap& normals, Arrangement& gm,
                        const NamedParameters& np = cp::default_values()) {
  using Face_normal_map = FaceNormalMap;
  using Graph_traits = boost::graph_traits<Graph>;
  using vertex_descriptor = typename Graph_traits::vertex_descriptor;
  using face_descriptor = typename Graph_traits::face_descriptor;

  auto points = CGAL::get(CGAL::vertex_point, g);
  std::unordered_map<vertex_descriptor, bool> vflags;
  auto flags = boost::make_assoc_property_map(vflags);
  Gaussian_map_builder<Graph, Face_normal_map, decltype(points), decltype(flags)> builder(gm, normals, points, flags);
  auto src = *CGAL::vertices(g).first;
  builder(g, src, true);
}

#endif
