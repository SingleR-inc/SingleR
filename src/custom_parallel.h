#ifndef CUSTOM_PARALLEL_H
#define CUSTOM_PARALLEL_H

#include <cmath>
#include <vector>
#include <thread>
#include <iostream>

extern int num_threads;

#define RATICATE_PARALLELIZE_UNKNOWN
#include "raticate/raticate.hpp"

template<class Function>
void parallelize(size_t n, Function function) {
    raticate::parallelize<double, int>(n, function, num_threads);
}

#define SINGLEPP_CUSTOM_PARALLEL parallelize

#endif
