cmake_minimum_required(VERSION 3.14...3.50)
set(project cgalex)
project(${project} VERSION "0.0.1" LANGUAGES CXX)

# Include windows toolchain path (if available)
if(CMAKE_TOOLCHAIN_FILE)
  include(${CMAKE_TOOLCHAIN_FILE})
endif(CMAKE_TOOLCHAIN_FILE)

# Define options
option(BUILD_SHARED_LIBS "Build shared libs instead of static libs" ON)
# option(${project}_BUILD_PYBINDINGS "Build python bindings" OFF)
option(CGALEX_USE_STATIC_LIBS "Link with static libraries" OFF)
# option(${project}_BUILD_DEV_DOC "Build the development documentation" OFF)
# option(${project}_BUILD_PYBINDINGS "Build python bindings" ON)
# option(${project}_BUILD_PYBINDINGS_DOC "Build the documentation of the Python bindings" ON)

# Options
set(${project}_PROJECT_AUTHORS "Efi Fogel")
set(${project}_PROJECT_EMAILS "efifogel@gmail.com")
set(${project}_PROJECT_URL "https://blah/foo")

# General flags

set(${project}_MODULES_REL_DIR cmake/modules)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/${${project}_MODULES_REL_DIR})
include(${project}Macros)
include(${project}Common)

# Specific options
option(${project}_TRACE "Trace" OFF)
option(${project}_WITH_VISUAL "With visual effects" OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
    STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()
set(default_build_type "Release")

# Binding options. Notice that The default value of option is OFF.
message(STATUS "VERSION=${CMAKE_PROJECT_VERSION}")
message(STATUS "SOVERSION=${CMAKE_PROJECT_VERSION_MAJOR}")

# Add src subdirectory
add_subdirectory(src)

include (InstallRequiredSystemLibraries)
#set (CPACK_RESOURCE_FILE_LICENSE
#     "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set(CPACK_GENERATOR TGZ)
set(CPACK_PACKAGE_VENDOR "Xanadu")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "1")
set(CPACK_PACKAGE_VERSION_PATCH "0")
include(CPack)
