#ifndef GAUSSIAN_MAP_3_H
#define GAUSSIAN_MAP_3_H

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Arr_default_overlay_traits.h>

#include "cgalex/Adder.h"

#include "arr_gaussian_map.h"
#include "Polyhedron_builder.h"

namespace pmp = CGAL::Polygon_mesh_processing;

template <typename Mesh, typename InputIterator>
Mesh minkowski_sum_3(InputIterator begin1, InputIterator end1, InputIterator begin2, InputIterator end2) {
  CGAL::Arr_face_overlay_traits<Arrangement, Arrangement, Arrangement, Adder> overlay_traits;
  Mesh ms;
  bool first = true;
  for (auto it1 = begin1; it1 != end1; ++it1) {
    for (auto it2 = begin2; it2 != end2; ++it2) {
      Arrangement gm;
      CGAL::overlay(*it1, *it2, gm, overlay_traits);
      using Halfedge_ds = typename Mesh::HalfedgeDS;
      Polyhedron_builder<Halfedge_ds, Arrangement> surface(gm);
      if (first) {
        ms.delegate(surface);
        pmp::triangulate_faces(ms);
        first = false;
        continue;
      }
      Mesh mesh;
      mesh.delegate(surface);
      pmp::triangulate_faces(mesh);
      pmp::corefine_and_compute_union(ms, mesh, ms);
    }
  }
  return ms;
}

#endif
