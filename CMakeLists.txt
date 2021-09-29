cmake_minimum_required(VERSION 3.20)
project(1_1)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(1_1
        1.1.cpp
        ParseLatexFunc.h
        SqMatrix.h
        TeX_convertible.h)
