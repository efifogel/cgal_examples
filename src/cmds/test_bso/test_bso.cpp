#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/container/small_vector.hpp>

#include <CGAL/Basic_viewer.h>
#include <CGAL/boost/graph/selection.h>
#include <CGAL/draw_surface_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/remesh_planar_patches.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Surface_mesh.h>

#include "CGAL/boolean_operations_3.h"

#ifdef CGAL_LINKED_WITH_TBB
using Concurrency_tag = CGAL::Parallel_tag;
#else
using Concurrency_tag = CGAL::Sequential_tag;
#endif

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;

using Mesh = CGAL::Surface_mesh<Kernel::Point_3>;
using vertex_descriptor = boost::graph_traits<Mesh>::vertex_descriptor;
using halfedge_descriptor = boost::graph_traits<Mesh>::halfedge_descriptor;
using edge_descriptor = boost::graph_traits<Mesh>::edge_descriptor;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;
using Triangle = boost::container::small_vector<std::size_t, 3>;

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

template <typename Mesh>
void apply(const Mesh& mesh1, const Mesh& mesh2, Mesh& result) {
  std::vector<Kernel::Point_3> points1, points2, points;
  std::vector<Triangle> triangles1, triangles2, triangles;
  PMP::polygon_mesh_to_polygon_soup(mesh1, points1, triangles1);
  PMP::polygon_mesh_to_polygon_soup(mesh2, points2, triangles2);
  CGAL::compute_difference<Concurrency_tag>(points1, triangles1, points2, triangles2, points, triangles);
  PMP::orient_polygon_soup(points, triangles);
  PMP::polygon_soup_to_polygon_mesh(points, triangles, result);
  PMP::stitch_borders(result, params::apply_per_connected_component(true));

  auto is_valid = result.is_valid();
  if (! is_valid) std::cerr << "The mesh is not valid\n";

  auto self_intersect = PMP::does_self_intersect(result);
  if (self_intersect) std::cerr << "The mesh self intersects\n";

  auto is_closed = CGAL::is_closed(result);
  if (! is_closed) std::cerr << "The mesh is not closed\n";
}

//! Main entry
int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cerr << "Usage: test_bso mesh1 mesh2 mesh3\n";
    return -1;
  }

  const auto* filename1 = argv[1];
  const auto* filename2 = argv[3];
  const auto* filename3 = argv[2];

  Mesh mesh1, mesh2, mesh3;
  if (! PMP::IO::read_polygon_mesh(filename1, mesh1)) {
    std::cerr << "Invalid input mesh 1" << std::endl;
    return 1;
  }
  if (! PMP::IO::read_polygon_mesh(filename2, mesh2)) {
    std::cerr << "Invalid input mesh 2" << std::endl;
    return 1;
  }
  if (! PMP::IO::read_polygon_mesh(filename3, mesh3)) {
    std::cerr << "Invalid input mesh 3" << std::endl;
    return 1;
  }

  CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso;
  gso.ignore_all_vertices(true);
  gso.ignore_all_edges(true);
  gso.colored_face = [](const Mesh&, face_descriptor) -> bool { return true; };
  gso.face_color = [] (const Mesh&, face_descriptor fh) -> CGAL::IO::Color {
    if (fh == boost::graph_traits<Mesh>::null_face()) return CGAL::IO::Color(100, 125, 200);
    return get_random_color(CGAL::get_default_random());
  };
  // CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso1;
  // gso1.ignore_all_vertices(true);
  // gso1.ignore_all_edges(true);
  // gso1.colored_face = [](const Mesh&, face_descriptor) -> bool { return true; };
  // gso1.face_color = [](const Mesh&, face_descriptor fh) -> CGAL::IO::Color { return CGAL::IO::Color(0, 0, 255); };
  // CGAL::Graphics_scene_options<Mesh, vertex_descriptor, edge_descriptor, face_descriptor> gso2;
  // gso2.ignore_all_vertices(true);
  // gso2.ignore_all_edges(true);
  // gso2.colored_face = [](const Mesh&, face_descriptor) -> bool { return true; };
  // gso2.face_color = [](const Mesh&, face_descriptor fh) -> CGAL::IO::Color { return CGAL::IO::Color(255, 0, 0); };
  // CGAL::Graphics_scene scene;
  // CGAL::add_to_graphics_scene(mesh1, scene, gso1);
  // CGAL::add_to_graphics_scene(mesh2, scene, gso2);
  // CGAL::draw_graphics_scene(scene);

  Mesh result1, result2;
  apply(mesh1, mesh2, result1);
  apply(result1, mesh3, result2);
  CGAL::draw(result2, gso, "result");

  // Mesh temp;
  // PMP::remesh_planar_patches(result2, temp);
  // std::swap(result2, temp);
  // CGAL::IO::write_polygon_mesh("result.off", result, CGAL::parameters::stream_precision(17));

  return 0;
}
