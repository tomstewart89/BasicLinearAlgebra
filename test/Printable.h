#pragma once

#include <cstddef>

class Print;

class Printable
{
  public:
    virtual size_t printTo(Print& p) const = 0;
};
