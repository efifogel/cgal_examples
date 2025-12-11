# Useful for macro names
string(TOUPPER ${project} project_upper)

# Use C++17 for this directory and its sub-directories.
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-fno-builtin)
endif()

# Set minimal required versions
set(${project}_PY_MIN_VERSION 3.8)
set(${project}_BOOST_MIN_VERSION 1.70.0)

################################################################################

# Find required libraries

if (${project}_USE_STATIC_LIBS)
  set(Boost_USE_STATIC_LIBS ON)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
  set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++ -static")
  set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
else()
  set(Boost_USE_STATIC_LIBS OFF)
endif()

# Find Boost
find_package(Boost REQUIRED COMPONENTS headers)

# Find Qt6 package
if (${project}_WITH_VISUAL)
  find_package(Qt6 QUIET COMPONENTS Core Gui Widgets)
  if (NOT Qt6_FOUND)
    message(WARNING "Qt6 not found; all Qt6-dependant binding will be excluded.")
    set(${project}_WITH_VISUAL OFF CACHE BOOL "" FORCE)
  endif()
endif()

# Find CGAL
if (${project}_WITH_VISUAL)
  find_package(CGAL REQUIRED COMPONENTS Core Qt6)
  if (NOT CGAL_Qt6_FOUND)
    message(WARNING "CGAL_Qt6 not found; all Qt6-dependant binding will be excluded.")
    set(${project}_WITH_VISUAL OFF CACHE BOOL "" FORCE)
  endif()
else()
  find_package(CGAL REQUIRED COMPONENTS Core)
endif()

# find_package(Eigen3 3.3.0 QUIET)
# include(CGAL_Eigen3_support)
# if (NOT TARGET CGAL::Eigen3_support)
#   message(FATAL_ERROR "This program requires the Eigen3 library, and will not be compiled.")
# endif()

################################################################################

# Add some compiler options
if (${project}_WITH_VISUAL)
  add_definitions(-D${project_upper}_HAS_VISUAL=)
endif()

if (MSVC)
  add_definitions(-D${project_upper}_HAS_DECLSPEC)
#   add_definitions(-GR)
#   add_definitions(-vmg)
endif()

# Use dynamic link for all ${project} libraries
add_definitions(-D${project_upper}_ALL_DYN_LINK)

if (${project}_TRACE)
  add_definitions(-D${project_upper}_TRACE)
endif()

# Installation
if (${project}_WIN32_CMAKE_ON_CYGWIN)
  exec_program(cygpath ARGS -w "${CMAKE_INSTALL_PREFIX}"
               OUTPUT_VARIABLE CMAKE_INSTALL_PREFIX2)
  file(TO_CMAKE_PATH ${CMAKE_INSTALL_PREFIX2} CMAKE_INSTALL_PREFIX)
endif()

# Offer the user the choice of overriding the installation directories
set (${project}_INSTALL_LIB_DIR lib CACHE PATH
  "Installation directory for libraries")
set (${project}_INSTALL_BIN_DIR bin CACHE PATH
  "Installation directory for executables")
set (${project}_INSTALL_INC_DIR include CACHE PATH
  "Installation directory for header files")
if(WIN32 AND NOT CYGWIN)
  set(${project}_DEF_INSTALL_CMAKE_DIR cmake)
else()
  set(${project}_DEF_INSTALL_CMAKE_DIR lib/cmake/${project})
endif()
set(${project}_INSTALL_CMAKE_DIR ${${project}_DEF_INSTALL_CMAKE_DIR} CACHE PATH
  "Installation directory for CMake files")

# Make relative paths absolute (needed later on)
foreach(p LIB BIN INC CMAKE PYTHON)
  set(var ${project}_INSTALL_${p}_DIR)
  if(NOT IS_ABSOLUTE "${${var}}")
    set(${var} "${CMAKE_INSTALL_PREFIX}/${${var}}")
  endif()
endforeach()

# Propagate
add_subdirectory(libs)
add_subdirectory(cmds)

include(CMakePackageConfigHelpers)

configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/${project}Config.cmake.in
    ${CMAKE_BINARY_DIR}/${project}Config.cmake
    INSTALL_DESTINATION lib/cmake/${project}
)

write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/${project}ConfigVersion.cmake"
  VERSION ${version}
  COMPATIBILITY AnyNewerVersion)

install(FILES
    "${CMAKE_BINARY_DIR}/${project}Config.cmake"
  DESTINATION lib/cmake/${project})
