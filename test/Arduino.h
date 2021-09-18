#pragma once

#include <iomanip>
#include <sstream>

struct Print
{
    std::stringstream buf;

    template <typename T>
    void print(const T &obj)
    {
        buf << obj;
    }

    void begin(int)
    {
        buf << std::fixed << std::showpoint << std::setprecision(2);
        buf.str("");
    }

} Serial;

using std::max;
