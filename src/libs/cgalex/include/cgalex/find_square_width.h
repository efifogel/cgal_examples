#ifndef CGALEX_FIND_SQUARE_WIDTH_H
#define CGALEX_FIND_SQUARE_WIDTH_H

#include "cgalex/square_distance.h"

template <typename Arrangement, typename Kernel>
std::pair<typename Kernel::FT, typename Kernel::Direction_3>
find_square_width(const Arrangement& gm, const Kernel& kernel) {
  typename Kernel::FT sw;
  typename Kernel::Direction_3 direction;
  auto it = gm.vertices_begin();
  for (; it != gm.vertices_end(); ++it) {
    if (it->degree() < 3) continue;
    sw = square_distance(it->incident_halfedges(), kernel);
    direction = it->point();
    break;
  }
  for (++it; it != gm.vertices_end(); ++it) {
    if (it->degree() < 3) continue;
    auto tmp = square_distance(it->incident_halfedges(), kernel);
    if (tmp < sw) {
      sw = tmp;
      direction = it->point();
    }
  }
  return std::make_pair(sw, direction);
}

#endif
