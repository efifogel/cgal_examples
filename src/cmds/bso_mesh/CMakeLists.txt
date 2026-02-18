set (component bso_mesh)
string(TOUPPER ${component} component_upper)

# Add source files
set(${component_upper}_SOURCE_FILES ${${component_upper}_SOURCE_FILES} ${component}.cpp)

if (BUILD_SHARED_LIBS)
  add_definitions(-DCGALEX_ALL_DYN_LINK)
endif()

add_executable(${component} ${${component_upper}_SOURCE_FILES})
set_target_properties(${component} PROPERTIES OUTPUT_NAME ${component})

target_include_directories(${component} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>)

# Find required Boost components
find_package(Boost CONFIG ${BOOST_MIN_VERSION} REQUIRED COMPONENTS system filesystem program_options)

target_link_libraries(${component} PRIVATE
  CGAL::CGAL
  Boost::program_options
  Boost::filesystem
  cgalex)

if (${project}_WITH_VISUAL)
  target_link_libraries(${component} PRIVATE CGAL::CGAL_Qt6 CGAL::CGAL_Basic_viewer)
endif ()

# Don't skip the full RPATH for the build tree
# set_property(TARGET player PROPERTY SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
# set_property(TARGET player PROPERTY BUILD_WITH_INSTALL_RPATH TRUE)

# the RPATH to be used when installing

# don't add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
if (NOT CGALEX_USE_STATIC_LIBS)
  set_property(TARGET ${component} PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
endif()

install(TARGETS ${component}
  EXPORT CGALEXTargets
  RUNTIME DESTINATION ${CGALEX_INSTALL_BIN_DIR}
  LIBRARY DESTINATION ${CGALEX_INSTALL_LIB_DIR}
  ARCHIVE DESTINATION ${CGALEX_INSTALL_LIB_DIR})
