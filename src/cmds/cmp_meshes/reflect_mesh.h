#ifndef REFLECT_MESH_PM
#define REFLECT_MESH_PM

// #include <algorithm>

#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Mesh, typename Kernel>
Mesh reflect_mesh(const Mesh& m, const Kernel& kernel) {
  using Transformation = CGAL::Aff_transformation_3<Kernel>;
  auto rm = m;
  Transformation reflect(CGAL::SCALING, -1.0);
  PMP::transform(reflect, rm);
  // std::transform(rm.points_begin(), rm.points_end(), rm.points_begin(),
  //                [](const Point& p) -> Point { return CGAL::ORIGIN - (p - CGAL::ORIGIN); });
  PMP::reverse_face_orientations(rm);
  // using Vector = typename Mesh::Plane_3;
  // std::transform(rm.planes_begin(), rm.planes_end(), rm.planes_begin(), [](const Vector& v) -> Vector { return -v; });
  return rm;
}

#endif
