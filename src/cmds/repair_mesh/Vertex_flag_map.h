#ifndef VERTEX_FLAG_MAP_MAP_H
#define VERTEX_FLAG_MAP_MAP_H

template <typename Polyhedron>
class Vertex_flag_map {
  using Vertex_handle = typename Polyhedron::Vertex_handle;

public:
  using key_type = Vertex_handle;
  using value_type = bool;
  using reference = bool;
  using category = boost::read_write_property_map_tag;

  Vertex_flag_map() {}

  friend bool get(const Vertex_flag_map&, Vertex_handle vd)
  { return vd->flag(); }

  friend void put(const Vertex_flag_map&, Vertex_handle vd, bool flag)
  { vd->set_flag(flag); }
};

#endif
