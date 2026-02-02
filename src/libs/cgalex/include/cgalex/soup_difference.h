#ifndef CGALEX_SOUP_DIFFERENCE_H
#define CGALEX_SOUP_DIFFERENCE_H

#include <vector>

#include <boost/container/small_vector.hpp>

#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include "CGAL/boolean_operations_3.h"

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

/*! computes the difference of two polygon meshes using polygon soup operations.
 */
template <typename Mesh>
void soup_difference(const Mesh& mesh1, const Mesh& mesh2, Mesh& result) {
  using Triangle = boost::container::small_vector<std::size_t, 3>;
  using Point_pm = typename boost::property_map<Mesh, boost::vertex_point_t>::type;
  using Point_3 = typename boost::property_traits<Point_pm>::value_type;
#ifdef CGAL_LINKED_WITH_TBB
using Concurrency_tag = CGAL::Parallel_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

  std::vector<Point_3> points1, points2;
  std::vector<Triangle> triangles1, triangles2;
  PMP::polygon_mesh_to_polygon_soup(mesh1, points1, triangles1);
  PMP::polygon_mesh_to_polygon_soup(mesh2, points2, triangles2);

  std::vector<Point_3> points;
  std::vector<Triangle> triangles;
  CGAL::compute_difference<Concurrency_tag>(points1, triangles1, points2, triangles2, points, triangles);
  PMP::orient_polygon_soup(points, triangles);

  result.clear();
  PMP::polygon_soup_to_polygon_mesh(points, triangles, result);
  PMP::stitch_borders(result, params::apply_per_connected_component(true));
}

#endif
