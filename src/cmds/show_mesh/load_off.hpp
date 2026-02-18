#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <tuple>
#include <cctype>
#include <iomanip>

// Simple structures
struct Vertex {
  double x, y, z;
  // optional color (0..1)
  bool hasColor = false;
  double r=0, g=0, b=0, a=1.0;
};

struct Face {
  std::vector<int> indices;
  // optional color
  bool hasColor = false;
  double r=0, g=0, b=0, a=1.0;
};

struct Mesh {
  std::vector<Vertex> verts;
  std::vector<Face> faces;
};

// trim helper
static inline std::string trim(const std::string &s) {
  size_t a = 0;
  while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
  size_t b = s.size();
  while (b > a && std::isspace((unsigned char)s[b-1])) --b;
  return s.substr(a, b-a);
}

// read next non-empty, non-comment line (comments start with '#')
bool getline_noncomment(std::istream &in, std::string &out) {
  while (std::getline(in, out)) {
    out = trim(out);
    if (out.empty()) continue;
    if (!out.empty() && out[0] == '#') continue;
    return true;
  }
  return false;
}

// parse whitespace separated tokens from a stringstream
template<typename T>
bool read_val(std::istringstream &iss, T &val) {
  iss >> val;
  return bool(iss);
}

// parse double token that might be integer color [0..255] or float [0..1]
double parse_color_token(const std::string &tok) {
  // If token contains a dot or 'e'/'E' treat as float; otherwise treat as int if >= 0 && <= 255
  bool hasDot = (tok.find('.') != std::string::npos) || (tok.find('e') != std::string::npos) || (tok.find('E') != std::string::npos);
  if (hasDot) {
    try { return std::stod(tok); }
    catch(...) { return 0.0; }
  } else {
    // integer-like: map 0..255 -> 0..1
    try {
      int v = std::stoi(tok);
      if (v < 0) v = 0;
      if (v > 255) v = 255;
      return v / 255.0;
    } catch(...) {
      // fallback to stod
      try { return std::stod(tok); } catch(...) { return 0.0; }
    }
  }
}

// Main OFF loader
bool loadOFF(const std::string &filename, Mesh &mesh, std::string *errOut = nullptr) {
  std::ifstream ifs(filename.c_str());
  if (!ifs) {
    if (errOut) *errOut = "Cannot open file: " + filename;
    return false;
  }

  std::string line;
  if (!getline_noncomment(ifs, line)) {
    if (errOut) *errOut = "Empty file or only comments.";
    return false;
  }

  std::string header = line;
  bool coff = false; // COFF indicates colors
  if (header == "OFF") {
    coff = false;
  } else if (header == "COFF") {
    coff = true;
  } else {
    // some OFF files include a leading "OFF" after optional whitespace or comment; allow lowercase as well
    std::string hs = header;
    for (auto &c: hs) c = (char)std::toupper((unsigned char)c);
    if (hs == "OFF") coff = false;
    else if (hs == "COFF") coff = true;
    else {
      // Possibly header missing and the first line has counts; try to detect by checking if header begins with a digit
      bool firstIsDigit = !hs.empty() && (std::isdigit((unsigned char)hs[0]) || hs[0]=='-' || hs[0]=='+');
      if (!firstIsDigit) {
        if (errOut) *errOut = "Not a valid OFF file (missing OFF/COFF header). Found: " + header;
        return false;
      } else {
        // Treat header line as counts line
        // fall through after setting line to header for counts parsing
      }
    }
  }

  // Read counts line (vertices faces edges)
  std::string countsLine;
  if (header == "OFF" || header == "COFF") {
    if (!getline_noncomment(ifs, countsLine)) {
      if (errOut) *errOut = "Missing counts line after OFF header.";
      return false;
    }
  } else {
    // header contained counts already (handled above)
    countsLine = header;
  }

  std::istringstream issCounts(countsLine);
  int nV=0, nF=0, nE=0;
  if (!(issCounts >> nV >> nF)) {
    if (errOut) *errOut = "Failed to parse vertex/face counts from: " + countsLine;
    return false;
  }
  issCounts >> nE; // optional

  if (nV < 0 || nF < 0) {
    if (errOut) *errOut = "Invalid counts in header.";
    return false;
  }

  mesh.verts.clear();
  mesh.faces.clear();
  mesh.verts.reserve(nV);
  mesh.faces.reserve(nF);

  // Read vertices
  for (int i = 0; i < nV; ++i) {
    if (!getline_noncomment(ifs, line)) {
      if (errOut) *errOut = "Unexpected EOF while reading vertices.";
      return false;
    }
    std::istringstream iss(line);
    Vertex v;
    if (!(iss >> v.x >> v.y >> v.z)) {
      if (errOut) *errOut = "Failed to read vertex coordinates on line: " + line;
      return false;
    }
    // optionally colors if COFF or extra tokens present
    std::string tok;
    if (coff) {
      // Expect r g b [a]
      if ((iss >> tok)) {
        v.hasColor = true;
        v.r = parse_color_token(tok);
        if (iss >> tok) v.g = parse_color_token(tok); else v.g = 0.0;
        if (iss >> tok) v.b = parse_color_token(tok); else v.b = 0.0;
        if (iss >> tok) v.a = parse_color_token(tok); else v.a = 1.0;
      }
    } else {
      // If non-COFF but there are >3 tokens, it might include colors; attempt to read them
      if (iss >> tok) {
        // There is an extra token; attempt to parse as color sequence
        v.hasColor = true;
        v.r = parse_color_token(tok);
        if (iss >> tok) v.g = parse_color_token(tok); else v.g = 0.0;
        if (iss >> tok) v.b = parse_color_token(tok); else v.b = 0.0;
        if (iss >> tok) v.a = parse_color_token(tok); else v.a = 1.0;
      }
    }
    mesh.verts.push_back(v);
  }

  // Read faces
  for (int i = 0; i < nF; ++i) {
    if (!getline_noncomment(ifs, line)) {
      if (errOut) *errOut = "Unexpected EOF while reading faces.";
      return false;
    }
    std::istringstream iss(line);
    int nv;
    if (!(iss >> nv)) {
      if (errOut) *errOut = "Failed to read face vertex count on line: " + line;
      return false;
    }
    if (nv < 0) {
      if (errOut) *errOut = "Invalid face vertex count: " + std::to_string(nv);
      return false;
    }
    Face f;
    f.indices.resize(nv);
    for (int k = 0; k < nv; ++k) {
      if (!(iss >> f.indices[k])) {
        if (errOut) *errOut = "Not enough indices for face on line: " + line;
        return false;
      }
      if (f.indices[k] < 0 || f.indices[k] >= nV) {
        if (errOut) *errOut = "Face index out of range: " + std::to_string(f.indices[k]);
        return false;
      }
    }
    // optional color for face (r g b [a]) or other tokens (e.g., material)
    std::string tok;
    if (iss >> tok) {
      // try parse as color tokens; if they are numbers we'll treat them as color
      f.hasColor = true;
      f.r = parse_color_token(tok);
      if (iss >> tok) f.g = parse_color_token(tok); else f.g = 0.0;
      if (iss >> tok) f.b = parse_color_token(tok); else f.b = 0.0;
      if (iss >> tok) f.a = parse_color_token(tok); else f.a = 1.0;
    }

    mesh.faces.push_back(std::move(f));
  }

  return true;
}

// Simple demo program
int main(int argc, char **argv) {
  if (argc < 2) {
    std::cerr << "Usage: off_reader <file.off>\n";
    return 1;
  }

  std::string err;
  Mesh mesh;
  if (!loadOFF(argv[1], mesh, &err)) {
    std::cerr << "Failed to load OFF: " << err << "\n";
    return 2;
  }

  std::cout << "Loaded OFF file: " << argv[1] << "\n";
  std::cout << "Vertices: " << mesh.verts.size() << ", Faces: " << mesh.faces.size() << "\n\n";

  // Print a few vertices
  for (size_t i = 0; i < mesh.verts.size() && i < 5; ++i) {
    const Vertex &v = mesh.verts[i];
    std::cout << "v["<<i<<"] = ("<<v.x<<","<<v.y<<","<<v.z<<")";
    if (v.hasColor) std::cout << " color=("<<v.r<<","<<v.g<<","<<v.b<<","<<v.a<<")";
    std::cout << "\n";
  }

  // Print a few faces
  for (size_t i = 0; i < mesh.faces.size() && i < 5; ++i) {
    const Face &f = mesh.faces[i];
    std::cout << "f["<<i<<"] =";
    for (int idx : f.indices) std::cout << " " << idx;
    if (f.hasColor) std::cout << " color=("<<f.r<<","<<f.g<<","<<f.b<<","<<f.a<<")";
    std::cout << "\n";
  }

  return 0;
}
