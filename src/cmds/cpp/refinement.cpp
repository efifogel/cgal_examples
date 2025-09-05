#include <type_traits>
#include <utility>

// "EqualityComparable"
template <typename T, typename = void>
struct is_equality_comparable : std::false_type {};

template <typename T>
struct is_equality_comparable<T,
  std::void_t<decltype(std::declval<T>() == std::declval<T>()),
              decltype(std::declval<T>() != std::declval<T>())>> : std::true_type
{};

// "TotallyOrdered" refines "EqualityComparable"
template <typename T, typename = void>
struct is_totally_ordered : std::false_type {};

template <typename T>
struct is_totally_ordered<T,
  std::enable_if_t<is_equality_comparable<T>::value,
    std::void_t<decltype(std::declval<T>() <  std::declval<T>()),
                decltype(std::declval<T>() <= std::declval<T>()),
                decltype(std::declval<T>() >  std::declval<T>()),
                decltype(std::declval<T>() >= std::declval<T>())>>> : std::true_type
{};

int main() {
  static_assert(is_equality_comparable<int>::value, "int should work");
  static_assert(is_totally_ordered<int>::value, "int is totally ordered");
  // static_assert(is_totally_ordered<void>::value); // error
}
