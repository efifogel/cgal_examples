#ifndef CGALEX_COREFINE_UNION_H
#define CGALEX_COREFINE_UNION_H

#include <vector>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>

#include "cgalex/soup_union.h"

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

/*! computes the union of two polygon meshes using corefine and compute union.
 */
template <typename Mesh>
void corefine_union(Mesh& mesh1, Mesh& mesh2, Mesh& result, std::size_t verbose_level = 0) {
  PMP::Corefinement::Non_manifold_output_visitor<Mesh> visitor(mesh1, mesh2);
  bool is_valid;
  try {
    is_valid = PMP::corefine_and_compute_union(mesh1, mesh2, result,
                                               params::visitor(visitor).throw_on_self_intersection(true));
  }
  catch (const PMP::Corefinement::Self_intersection_exception& e) {
    if (verbose_level > 0) {
      std::cout << "UNION: Caught an exception: " << e.what() << std::endl;
    }

    soup_union(mesh1, mesh2, result);
    return;
  }

  if (is_valid) {
    auto num_vertices_created = PMP::duplicate_non_manifold_vertices(result);
    if (verbose_level > 1) {
      if (num_vertices_created > 0) {
        std::cout << "UNION: Created " << num_vertices_created << " to resolve non-manifold vertices" << std::endl;
      }
      auto self_intersect = PMP::does_self_intersect(result);
      if (self_intersect)
        std::cout << "UNION: The output self intersect as a result of duplicating non-manifold vertices" << std::endl;
    }
  }
  else {
    if (verbose_level > 0) {
      std::cout << "UNION: Failed - resorting to handling a non-manifold result\n";
    }
    using Point_pm = typename boost::property_map<Mesh, boost::vertex_point_t>::type;
    using Point_3 = typename boost::property_traits<Point_pm>::value_type;
    std::vector<Point_3> points;
    std::vector<std::array<std::size_t, 3>> polygons;
    visitor.extract_union(points, polygons);
    PMP::orient_polygon_soup(points, polygons);
    PMP::polygon_soup_to_polygon_mesh(points, polygons, result);
    PMP::stitch_borders(result, params::apply_per_connected_component(true));
  }
}

#endif
