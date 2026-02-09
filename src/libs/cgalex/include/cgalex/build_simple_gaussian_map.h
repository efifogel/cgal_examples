#ifndef CGALEX_BUILD_SIMPLE_GAUSSIAN_MAP_H
#define CGALEX_BUILD_SIMPLE_GAUSSIAN_MAP_H

#include "cgalex/Simple_gaussian_map_builder.h"

template <typename Graph, typename FaceNormalMap>
void build_simple_gaussian_map(const Graph& g, const FaceNormalMap& normals,
                               Arrangement& gm) {
  std::cout << "before\n";
  auto points = CGAL::get(CGAL::vertex_point, g);
  Simple_gaussian_map_builder<Graph, FaceNormalMap, decltype(points)> builder(gm, normals, points);
  auto src = *CGAL::vertices(g).first;
  builder(g, src, true);
}

#endif
