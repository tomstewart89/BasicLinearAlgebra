#pragma once

#include <iostream>

struct Print {
  template <typename T> void print(const T &obj) { std::cout << obj; }
} Serial;
