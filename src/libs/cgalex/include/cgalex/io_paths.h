#ifndef CGALEX_IO_PATHS_HPP
#define CGALEX_IO_PATHS_HPP

#include <boost/filesystem.hpp>

#include "cgalex/Paths.h"

namespace fi = boost::filesystem;

namespace std {

/*! Export `fi::path` to an output stream.
 */
template <typename OutputStream>
OutputStream& operator<<(OutputStream& os, const Path& path) {
  os << path.string();
  return os;
}

/*! Export Paths to an output stream.
 */
template <typename OutputStream>
OutputStream& operator<<(OutputStream& os, const Paths& paths) {
  for (const auto& path : paths) {
    os << path.string();
    if (path != paths.back()) os << " ";
  }
  return os;
}

}

#endif
