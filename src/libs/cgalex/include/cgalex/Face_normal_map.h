#ifndef CGALEX_FACE_NORMAL_MAP_MAP_H
#define CGALEX_FACE_NORMAL_MAP_MAP_H

template <typename Polyhedron>
class Face_normal_map {
  using Plane_3 = typename Polyhedron::Plane_3;
  using Facet_handle = typename Polyhedron::Facet_handle;

public:
  using key_type = Facet_handle;
  using value_type = Plane_3;
  using reference = Plane_3;
  using category = boost::read_write_property_map_tag;

  Face_normal_map() {}

  friend value_type get(const Face_normal_map&, Facet_handle f)
  { return f->plane(); }

  friend void put(const Face_normal_map&, key_type f, const value_type& val)
  { f->plane() = val; }
};

#endif
