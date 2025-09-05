#include <type_traits>

template <typename T>
typename std::enable_if<std::is_integral<T>::value, T>::type
add_one(T x) {
  return x + 1;
}

int main() {
  int a = add_one(5);          // OK
  // double d = add_one(3.14); // Error: not integral
  return 0;
}
