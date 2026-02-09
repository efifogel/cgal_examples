#ifndef CGALEX_GAUSSIAN_MAP_H
#define CGALEX_GAUSSIAN_MAP_H

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include "cgalex/merge_coplanar_faces.h"
#include "cgalex/build_gaussian_map.h"
#include "cgalex/Face_normal_map.h"
#include "cgalex/Vertex_flag_map.h"

template <typename GaussianMap, typename Mesh, typename Kernel>
GaussianMap gaussian_map(Mesh& mesh, const Kernel& kernel) {
  GaussianMap gm;
  auto np = CGAL::parameters::geom_traits(kernel);
  Face_normal_map<Mesh> normals;
  CGAL::Polygon_mesh_processing::compute_face_normals(mesh, normals, np);
  merge_coplanar_faces(mesh, normals, np);
  Vertex_flag_map<Mesh> vertex_flags;
  auto npv = CGAL::parameters::vertex_index_map(vertex_flags);
  build_gaussian_map(mesh, normals, gm, npv);
  return gm;
}

#endif
