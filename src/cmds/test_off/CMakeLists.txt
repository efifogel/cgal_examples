cmake_minimum_required(VERSION 3.1...3.23)
project(test_off)
find_package(Qt6 QUIET COMPONENTS Widgets)
find_package(CGAL QUIET COMPONENTS Core OPTIONAL_COMPONENTS Qt6)
if (NOT CGAL_FOUND)
  message(FATAL_ERROR "This program requires the CGAL library, and will not be compiled.")
endif()

include_directories(
  "../chapter_2"
  "../chapter_3"
)

create_single_source_cgal_program("test_off.cpp")

if(CGAL_Qt6_FOUND)
  target_link_libraries(test_off PUBLIC CGAL::CGAL_Basic_viewer)
else()
  message(STATUS "NOTICE: The drawings in some examples are is disabled, as they depends on Qt, and Qt is not supported.")
endif()
