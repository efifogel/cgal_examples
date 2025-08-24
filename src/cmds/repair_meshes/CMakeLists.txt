cmake_minimum_required(VERSION 3.1...3.23)
project(repair_meshes LANGUAGES CXX)
# Find Boost with program_options
find_package(Boost REQUIRED COMPONENTS program_options)
find_package(Qt6 QUIET COMPONENTS Widgets)
find_package(CGAL QUIET COMPONENTS Core OPTIONAL_COMPONENTS Qt6)
if (NOT CGAL_FOUND)
  message(FATAL_ERROR "This program requires the CGAL library, and will not be compiled.")
endif ()
if (NOT CGAL_Qt6_FOUND)
  message(FATAL_ERROR "This program requires the CGAL_Qt6 library, and will not be compiled.")
endif()

find_package(Eigen3 3.2.0 QUIET) #(requires 3.2.0 or greater)
include(CGAL_Eigen3_support)
if (NOT TARGET CGAL::Eigen3_support)
  message(FATAL_ERROR "This program requires the Eigen3 library, and will not be compiled.")
endif()

# create_single_source_cgal_program("repair_meshes.cpp")
add_executable(repair_meshes repair_meshes.cpp)

target_link_libraries(repair_meshes PUBLIC CGAL::CGAL_Basic_viewer CGAL::Eigen3_support Boost::program_options)
