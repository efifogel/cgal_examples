#ifndef CGALEX_GAUSSIAN_MAP_SM_H
#define CGALEX_GAUSSIAN_MAP_SM_H

#include <boost/graph/graph_traits.hpp>

#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include "cgalex/merge_coplanar_faces.h"
#include "cgalex/build_simple_gaussian_map.h"

template <typename GaussianMap, typename Mesh, typename Kernel>
GaussianMap gaussian_map_sm(Mesh& mesh, const Kernel& kernel) {
  using Vector_3 = typename Kernel::Vector_3;
  using Fd = typename boost::graph_traits<Mesh>::face_descriptor;
  auto np = CGAL::parameters::geom_traits(kernel);
  auto normals = mesh.template add_property_map<Fd, Vector_3>("f:normals", CGAL::NULL_VECTOR).first;
  CGAL::Polygon_mesh_processing::compute_face_normals(mesh, normals, np);
  merge_coplanar_faces(mesh, normals, np);
  GaussianMap gm;
  build_simple_gaussian_map(mesh, normals, gm);
  return gm;
}

#endif
