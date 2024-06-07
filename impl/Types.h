#pragma once

namespace BLA
{
// This namespace exists because the header "typetraits" is not implemented in every Arduino environment.
namespace Types
{
template<class T, class U> struct is_same       { static constexpr bool value = false; };
template<class T>          struct is_same<T, T> { static constexpr bool value = true;  };

template<class T> struct remove_const { typedef T type; };
template<class T> struct remove_const<const T> { typedef T type; };

template<class T>
struct is_floating_point
{
    static constexpr bool value =
        is_same<float,       typename remove_const<T>::type>::value ||
        is_same<double,      typename remove_const<T>::type>::value ||
        is_same<long double, typename remove_const<T>::type>::value;
};

template<class T>
struct is_signed_integer
{
    static constexpr bool value =
        is_same<signed char,      typename remove_const<T>::type>::value ||
        is_same<signed short,     typename remove_const<T>::type>::value ||
        is_same<signed int,       typename remove_const<T>::type>::value ||
        is_same<signed long,      typename remove_const<T>::type>::value ||
        is_same<signed long long, typename remove_const<T>::type>::value;
};

template<class T>
struct is_signed
{
    static constexpr bool value =
        is_floating_point<T>::value ||
        is_signed_integer<T>::value;
};

template<bool, typename T = void> struct enable_if {};
template<typename T> struct enable_if<true, T> { typedef T type; };
}
} // namespace BLA

