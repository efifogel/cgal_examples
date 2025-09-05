#include <type_traits>
#include <vector>
#include <list>

template <typename C>
// typename std::enable_if<
//   std::is_same<C, std::vector<typename C::value_type>>::value ||
//   std::is_same<C, std::list<typename C::value_type>>::value
// >::type
typename std::enable_if<
  std::is_same<C, std::vector<typename C::value_type>>::value ||
  std::is_same<C, std::list<typename C::value_type>>::value, void
>::type
process_container(const C& c) {
  // Works for vector and list only
}

int main() {
  std::vector<int> v = {1,2,3};
  process_container(v); // OK

  // std::set<int> s = {1,2,3};
  // process_container(s); // Error
}
