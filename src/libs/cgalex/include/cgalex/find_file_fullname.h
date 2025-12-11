#ifndef CGALEX_FIND_FILE_FULLNAME_HPP
#define CGALEX_FIND_FILE_FULLNAME_HPP

#include <string>

#include <boost/filesystem.hpp>

namespace fi = boost::filesystem;

// CGALEX_BEGIN_NAMESPACE

//! \brief finds the input file.
template <typename Paths>
inline std::string find_file_fullname(const Paths& paths, const std::string& const_filename) {
  std::string filename(const_filename);

  // CGALEX_assertion(!filename.empty());

#if (defined _MSC_VER)
  // Convert the ROOT from cygwin path to windows path, if relevant:
  auto cygdrive = filename.substr(0, 10);
  if (cygdrive == std::string("/cygdrive/")) {
    filename.erase(0, 10);
    filename.insert(1, ":");
  }
#endif

  fi::path file_path(filename);
  if (file_path.is_absolute()) {
    if (fi::exists(file_path)) return file_path.string();
    return std::string();
  }

  for (const auto& path : paths) {
    auto full_file_path = path / file_path;
    if (fi::exists(full_file_path) && fi::is_regular_file(full_file_path)) return full_file_path.string();
  }
  return std::string();
}

// CGALEX_END_NAMESPACE

#endif
