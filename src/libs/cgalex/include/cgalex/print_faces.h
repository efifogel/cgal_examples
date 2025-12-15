#ifndef CGALEX_PRINT_FACES
#define CGALEX_PRINT_FACES

#include <boost/graph/graph_traits.hpp>

template <typename Mesh>
void print_face(typename boost::graph_traits<Mesh>::face_descriptor f, const Mesh& mesh) {
  std::cout << "    Switch {\n"
    "      whichChoice 0\n"
    "      children Shape {\n"
    "        appearance Appearance {\n"
    "          material Material { }\n"
    "        }\n"
    "        geometry IndexedLineSet {\n"
    "          coord Coordinate {\n"
    "            point [\n";

  std::size_t ci = 0;
  double color[3] = {1.0, 1.0, 1.0};

  // Get a halfedge bounding the face
  auto h = mesh.halfedge(f);

  // Iterate CCW around the face
  auto h_start = h;
  std::size_t cnt = 0;
  do {
    auto v = mesh.target(h);
    const auto& p = mesh.point(v);
    std::cout << "                   " << p.x() << " " << p.y() << " " << p.z() << std::endl;
    h = mesh.next(h);
    ++cnt;
  } while (h != h_start);

  std::cout << "                  ]\n"
    "          }\n"
    "          coordIndex [\n";
  std::cout << "                      ";
  for (auto i = 0; i < cnt; ++i) std::cout << i << " ";
  std::cout << "0 -1\n";
  std::cout << "                     ]\n"
    "          colorPerVertex FALSE\n"
    "          color Color {\n"
    "            color[\n";
  std::cout << "                  " << color[0] << " " << color[1] << " " << color[2] << "\n";
  std::cout << "                 ]\n"
    "          }\n"
    "        }\n"
    "      }\n"
    "    }\n";
  ci = (ci + 1) % 3;
  color[ci] -= 0.33;
  if (color[ci] < 0.0) color[ci] = 1.0;
  if ((color[0] < 0.1) && (color[1] < 0.1) && (color[2] < 0.1))
    color[0] = color[1] = color[2] = 1.0;
}

template <typename Mesh, typename Map, typename Direction>
void print_faces(const Mesh& mesh, const Map& normals) {
  std::cout << "#VRML V2.0 utf8\n"
    "NavigationInfo { type [ \"EXAMINE\" \"ANY\" ] }\n"
    "Group {\n"
    "  children [\n";

  using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
  Kernel kernel;
  auto eq = kernel.equal_3_object();

  for (auto f : mesh.faces()) {
    const auto& normal = get(normals, f);
    std::cout << "    # Face: " << f << ", " << normal << "\n";
    print_face(f, mesh);
  }
  std::cout << "   ]\n"
    "}  \n";
}

#endif
