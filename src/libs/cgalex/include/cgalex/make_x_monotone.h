#ifndef CGALEX_MAKE_X_MONOTONE_H
#define CGALEX_MAKE_X_MONOTONE_H

template <typename Vector_3, typename Traits, typename OutputIterator>
OutputIterator make_x_monotone(const Vector_3& n1, const Vector_3& n2,
                               OutputIterator oi, const Traits& traits) {
  auto ctr_point = traits.construct_point_2_object();
  auto ctr_cv = traits.construct_curve_2_object();
  auto cv = ctr_cv(ctr_point(n1.direction()), ctr_point(n2.direction()));
  *oi++ = traits.make_x_monotone_2_object()(cv, oi);
  return oi;
}

#endif
