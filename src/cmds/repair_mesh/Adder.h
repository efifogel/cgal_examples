#ifndef ADDER_H
#define ADDER_H

#include <CGAL/Kernel_traits.h>

struct Adder {
  template <typename Point>
  Point operator()(const Point& p, const Point& q) const {
    using Kernel = typename CGAL::Kernel_traits<Point>::Kernel;
    using Vector_3 = typename Kernel::Vector_3;
    return q + Vector_3(CGAL::ORIGIN, p);
  }
};

#endif
