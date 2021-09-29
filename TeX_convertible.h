#ifndef TEX_CONVERTIBLE_H
#define TEX_CONVERTIBLE_H
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <limits>
using namespace std;

class TeX_convertible{
    protected:
        virtual std::string convert() const = 0;
};
#endif


