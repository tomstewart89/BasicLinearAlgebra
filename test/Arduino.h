#pragma once

#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <sstream>

#include "Printable.h"

struct Print
{
    std::stringstream buf;

    template <typename T>
    typename std::enable_if<!std::is_base_of<Printable, T>::value, size_t>::type
    print(const T& obj)
    {
        buf << obj;
        return 0;
    }

    template <typename T>
    typename std::enable_if<!std::is_base_of<Printable, T>::value, size_t>::type
    println(const T& obj)
    {
        buf << obj << std::endl;
        return 0;
    }

    size_t print(const Printable& x)
    {
        x.printTo(*this);
        return 0;
    }

    size_t println(const Printable& x)
    {
        x.printTo(*this);
        println("");
        return 0;
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
