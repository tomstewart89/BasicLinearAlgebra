#pragma once

#include <algorithm>
#include <iomanip>
#include <sstream>

struct Print
{
    std::stringstream buf;

    template <typename T>
    void print(const T& obj)
    {
        buf << obj;
    }

    template <typename T>
    void println(const T& obj)
    {
        buf << obj << std::endl;
    }

    void begin(int)
    {
        buf << std::fixed << std::showpoint << std::setprecision(2);
        buf.str("");
    }

    Print& operator<<(std::ostream& (*pf)(std::ostream&))
    {
        buf << pf;
        return *this;
    }

} Serial;

using std::endl;
using std::max;
