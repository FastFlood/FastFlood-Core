#ifndef OSS_RASTERREDUCE_H
#define OSS_RASTERREDUCE_H

#include "map.h"
#include <atomic>
#include <mutex>
#include <algorithm>
#include "quickpool.h"

/**
 * @brief Generic map reduction function using a threadpool.
 *
 * @tparam T Data type of the map.
 * @tparam Callable Binary reduction operator.
 * @tparam CallableAtomic Atomic update function for thread-safe aggregation.
 * @param map Input map to reduce.
 * @param f The reduction function (e.g., std::plus).
 * @param g The atomic update function.
 * @param count Output parameter for the number of valid cells processed.
 * @return The reduced value.
 */
template <typename T, typename Callable, typename CallableAtomic>
inline double MapReduceByOperatorMultiThreadedCounted(const cTMap<T> *map, Callable && f, CallableAtomic && g, size_t & count)
{
    std::atomic<double> result = 0.0;
    std::atomic<size_t> count_t = 0;

    GetTPool().parallel_for(0, map->nrRows(), [&](size_t j) {
        double resultj = 0.0; // Note: Needs proper identity initialization
        size_t r = j;
        size_t counth = 0;
        for(size_t c = 0; c < map->nrCols(); c++) {
            if(!pcr::isMV(map->data[r][c])) {
                resultj = f(resultj, map->data[r][c]);
                counth++;
            }
        }
        count_t += counth;
        g(result, resultj);
    });

    GetTPool().waitext();
    count = count_t;
    return result.load();
}

/**
 * @brief Reduces a map to its minimum value.
 */
inline double MapMinimumRed(const cTMapf * map)
{
    size_t count = 0;
    std::mutex l;
    double min_val = 1e31;
    double min = MapReduceByOperatorMultiThreadedCounted(map,
        [&](const float &x, const float &y){return std::min(x,y);},
        [&](std::atomic<double>&x, const double &y){ l.lock(); x = std::min(x.load(), y); l.unlock();},
        count
    );
    if(count == 0) pcr::setMV((float&)min);
    return min;
}

/**
 * @brief Reduces a map to its maximum value.
 */
inline double MapMaximumRed(const cTMapf * map)
{
    size_t count = 0;
    std::mutex l;
    double max_val = -1e31;
    double max = MapReduceByOperatorMultiThreadedCounted(map,
        [&](const float &x, const float &y){return std::max(x,y);},
        [&](std::atomic<double>&x, const double &y){ l.lock(); x = std::max(x.load(), y); l.unlock();},
        count
    );
    if(count == 0) pcr::setMV((float&)max);
    return max;
}

/**
 * @brief Reduces a map to its average value.
 */
inline double MapAverageRed(const cTMapf * map)
{
    size_t count = 0;
    double total = MapReduceByOperatorMultiThreadedCounted(map,
        [](const float &x, const float &y){return x+y;},
        [](std::atomic<double>&x, const double &y){ x += y;},
        count
    );
    return total / std::max((size_t)1, count);
}

#endif
