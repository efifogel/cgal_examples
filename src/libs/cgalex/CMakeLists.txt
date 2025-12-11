set(component cgalex)
string(TOUPPER ${component} component_upper)

################################################################################

# Define directories

set(${component}_BINARY_SRC_DIR "${CMAKE_BINARY_DIR}/src")
set(${component}_BINARY_LIB_DIR "${${component}_BINARY_SRC_DIR}/libs/${component}")

# Define generated directories

set (${component}_GENERATED_BASE_DIR "${${component}_BINARY_LIB_DIR}")
set (${component}_GENERATED_DIR "${${component}_GENERATED_BASE_DIR}/generated")
set (${component}_GENERATED_LIB_DIR "${${component}_GENERATED_BASE_DIR}/lib/generated")
set (${component}_GENERATED_INC_DIR "${${component}_GENERATED_BASE_DIR}/include/generated")

file (MAKE_DIRECTORY ${${component}_GENERATED_DIR})

################################################################################

# Generate version file

configure_file(${project}_version.hpp.in ${${project}_GENERATED_INC_DIR}/${project}/version.hpp)

################################################################################

# Create custom commands

################################################################################

# Define source files

# Private headers

# Public headers
set(${component}_PUB_HDR_FILES ${${component}_PUB_HDR_FILES} include/${project}/find_file_fullname.h)
set(${component}_PUB_HDR_FILES ${${component}_PUB_HDR_FILES} include/${project}/io_paths.h)
set(${component}_PUB_HDR_FILES ${${component}_PUB_HDR_FILES} include/${project}/merge_coplanar_faces.h)
set(${component}_PUB_HDR_FILES ${${component}_PUB_HDR_FILES} include/${project}/Paths.h)
set(${component}_PUB_HDR_FILES ${${component}_PUB_HDR_FILES} include/${project}/retriangulate_faces.h)
set(${component}_PUB_HDR_FILES ${${component}_PUB_HDR_FILES} include/${project}/triangulate_faces.h)

source_group(src FILES ${${component}_SRC_FILES} ${${component}_PUB_HDR_FILES})

# Add the generated files

###############################################################################

# Add library

add_library(${component} INTERFACE)

target_sources(${component}
  PUBLIC
    FILE_SET pub_HEADERS # Name of the set, must be unique for a target
    TYPE HEADERS # This set describes headers
    BASE_DIRS include
    FILES ${${component}_PUB_HDR_FILES})

set_target_properties(${component} PROPERTIES
  VERSION "${CMAKE_PROJECT_VERSION}"
  SOVERSION "${CMAKE_PROJECT_VERSION_MAJOR}")

###############################################################################

# Find packages

###############################################################################

# Define header directories

target_include_directories(${component} INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)

###############################################################################

# Add defines

add_definitions(-D${project_upper}_${component_upper}_SOURCE)
add_definitions(-D${component_upper}_LIB)

# if (MSVC)
#   add_definitions(/wd4251)
# else()
#   add_definitions(-DGL_GLEXT_PROTOTYPES)
# endif()

################################################################################

# Add link options

################################################################################

# Install

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${component}Config.cmake.in
    ${CMAKE_BINARY_DIR}/${project}${component}Config.cmake
    INSTALL_DESTINATION lib/cmake/${project})

write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/${project}${component}ConfigVersion.cmake"
  VERSION ${version}
  COMPATIBILITY AnyNewerVersion)

install(TARGETS ${component} EXPORT ${component}Targets
  COMPONENT ${component}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
  RUNTIME DESTINATION bin
  FILE_SET pub_HEADERS DESTINATION include
  INCLUDES DESTINATION include)

install(EXPORT ${component}Targets
  FILE "${project}${component}Targets.cmake"
  NAMESPACE ${project}::
  DESTINATION lib/cmake/${project}
  COMPONENT ${component})

install(FILES
    "${CMAKE_BINARY_DIR}/${project}${component}Config.cmake"
    "${CMAKE_BINARY_DIR}/${project}${component}ConfigVersion.cmake"
  DESTINATION lib/cmake/${project}
  COMPONENT ${component})
