#ifndef HCI_ARRAY_MATH_H_
#define HCI_ARRAY_MATH_H_

#include <array>
#include <cstddef>
#include <iostream>

template <class T, std::size_t N>
std::array<T, N> operator+(
    const std::array<T, N> &lhs, const std::array<T, N> &rhs) {
  std::array<T, N> res;
  for (std::size_t i = 0; i < N; i++) res[i] = lhs[i] + rhs[i];
  return res;
}

template <class T, std::size_t N>
std::array<T, N> operator-(
    const std::array<T, N> &lhs, const std::array<T, N> &rhs) {
  std::array<T, N> res;
  for (std::size_t i = 0; i < N; i++) res[i] = lhs[i] - rhs[i];
  return res;
}

template <class T, std::size_t N>
std::array<double, N> operator*(const std::array<T, N> &lhs, const double rhs) {
  std::array<double, N> res;
  for (std::size_t i = 0; i < N; i++) res[i] = lhs[i] * rhs;
  return res;
}

template <class T, std::size_t N>
std::array<T, N> &operator+=(
    std::array<T, N> &lhs, const std::array<T, N> &rhs) {
  for (std::size_t i = 0; i < N; i++) lhs[i] += rhs[i];
  return lhs;
}

template <class T, std::size_t N>
std::array<T, N> &operator-=(
    std::array<T, N> &lhs, const std::array<T, N> &rhs) {
  for (std::size_t i = 0; i < N; i++) lhs[i] -= rhs[i];
  return lhs;
}

template <class T, std::size_t N>
std::array<T, N> square(const std::array<T, N> &arr) {
  std::array<T, N> res;
  for (std::size_t i = 0; i < N; i++) res[i] = arr[i] * arr[i];
  return res;
}

template <class T, std::size_t N>
T sum(const std::array<T, N> &arr) {
  T res = arr[0];
  for (std::size_t i = 1; i < N; i++) res += arr[i];
  return res;
}

template <class T, std::size_t N>
double norm(const std::array<T, N> &arr) {
  double res = arr[0] * arr[0];
  for (std::size_t i = 1; i < N; i++) res += arr[i] * arr[i];
  return sqrt(res);
}

template <class T, std::size_t N>
std::array<int, N> round(const std::array<T, N> &arr) {
  std::array<int, N> res;
  for (std::size_t i = 0; i < N; i++) {
    res[i] = static_cast<int>(ceil(arr[i] - 0.5));
  }
  return res;
}

template <class T, std::size_t N>
bool operator==(const std::array<T, N> &lhs, const T &rhs) {
  for (std::size_t i = 0; i < N; i++) {
    if (lhs[i] != rhs) return false;
  }
  return true;
}

template <class T, std::size_t N>
bool operator!=(const std::array<T, N> &lhs, const T &rhs) {
  return !(lhs == rhs);
}

template <class T, std::size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T, N> &rhs) {
  for (std::size_t i = 0; i < N - 1; i++) os << rhs[i] << ' ';
  os << rhs[N - 1];
  return os;
}

#endif