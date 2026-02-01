#ifndef GAUSSIAN_MAPS_H
#define GAUSSIAN_MAPS_H

#include <CGAL/convex_hull_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>
#include <CGAL/convex_decomposition_3.h>

#include "gaussian_map.h"

template <typename Mesh, typename OutputIterator, typename Kernel>
OutputIterator
gaussian_maps(Mesh& mesh, OutputIterator oi, const Kernel& kernel) {
  if (CGAL::is_strongly_convex_3(mesh, kernel)) {
    *oi++ = gaussian_map<Arrangement>(mesh, kernel);
    return oi;
  }
  CGAL::Nef_polyhedron_3<Kernel, CGAL::SNC_indexed_items> nef_mesh(mesh);
  CGAL::convex_decomposition_3(nef_mesh);
  auto ci = ++nef_mesh.volumes_begin();
  for (; ci != nef_mesh.volumes_end(); ++ci) {
    if (! ci->mark()) continue;
    Mesh mesh;
    nef_mesh.convert_inner_shell_to_polyhedron(ci->shells_begin(), mesh);
    *oi++ = gaussian_map<Arrangement>(mesh, kernel);
  }
  return oi;
}

#endif
