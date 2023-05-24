#ifndef SINGLEPP_MACROS_HPP
#define SINGLEPP_MACROS_HPP

/**
 * @file macros.hpp
 *
 * @brief Set common macros used throughout **singlepp**.
 *
 * @details
 * The `SINGLEPP_CUSTOM_PARALLEL` macro can be set to a function that specifies a custom parallelization scheme.
 * This function should be a template that accept three arguments:
 *
 * - `fun`, a lambda that accepts three arguments: 
 *   - `thread`, the thread number.
 *   - `start`, the index of the first job for this thread.
 *   - `length`, the number of jobs for this thread.
 * - `njobs`, an integer specifying the number of jobs.
 * - `nthreads`, an integer specifying the number of threads to use.
 *
 * The function should split `[0, njobs)` into any number of contiguous, non-overlapping intervals, and call `fun` on each interval, possibly in different threads.
 * The details of the splitting and evaluation are left to the discretion of the developer defining the macro. 
 * The function should only return once all evaluations of `fun` are complete.
 *
 * If `SINGLEPP_CUSTOM_PARALLEL` is set, the following macros are also set (if they are not already defined):
 *
 * - `TATAMI_CUSTOM_PARALLEL`, from the [**tatami**](https://ltla.github.io/tatami) library.
 *
 * This ensures that any custom parallelization scheme is propagated to all of **singlepp**'s dependencies.
 * If these libraries are used outside of **singlepp**, some care is required to ensure that the macros are consistently defined through the client library/application;
 * otherwise, developers may observe ODR compilation errors. 
 */

// Synchronizing all parallelization schemes.
#ifdef SINGLEPP_CUSTOM_PARALLEL

#ifndef TATAMI_CUSTOM_PARALLEL
#define TATAMI_CUSTOM_PARALLEL SINGLEPP_CUSTOM_PARALLEL
#endif

#endif

#endif
