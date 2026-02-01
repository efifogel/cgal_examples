#ifndef ARR_GAUSSIAN_MAP_H
#define ARR_GAUSSIAN_MAP_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arr_geodesic_arc_on_sphere_traits_2.h>
#include <CGAL/Arr_spherical_topology_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>

using Kernel = CGAL::Exact_predicates_exact_constructions_kernel;
using Point_3 = Kernel::Point_3;

using Geom_traits = CGAL::Arr_geodesic_arc_on_sphere_traits_2<Kernel>;
using Point = Geom_traits::Point_2;
using X_monotone_curve = Geom_traits::X_monotone_curve_2;
using Dcel = CGAL::Arr_face_extended_dcel<Geom_traits, Point_3>;
using Topol_traits = CGAL::Arr_spherical_topology_traits_2<Geom_traits,Dcel>;
using Arrangement = CGAL::Arrangement_on_surface_2<Geom_traits,Topol_traits>;
using Vertex_handle = Arrangement::Vertex_handle;
using Halfedge_handle = Arrangement::Halfedge_handle;
using Face_handle = Arrangement::Face_handle;

#endif
