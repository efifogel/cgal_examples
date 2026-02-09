#ifndef CGALEX_RETRIANGULE_COPLANAR_FACES_H
#define CGALEX_RETRIANGULE_COPLANAR_FACES_H

#include <CGAL/boost/graph/named_params_helper.h>

#include "cgalex/merge_coplanar_faces.h"
#include "cgalex/triangulate_faces.h"

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

/*! Retriangulate a mesh
 */
template <typename PolygonMesh, typename Map, typename NamedParameters = CGAL::parameters::Default_named_parameters>
void retriangulate_faces(PolygonMesh& mesh, const Map& normals, const NamedParameters& np = params::default_values()) {
  merge_coplanar_faces(mesh, normals, np);
  triangulate_faces(mesh, normals, np);
}

#endif
