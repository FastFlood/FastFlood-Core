#ifndef DATETIME_H
#define DATETIME_H

#include <iostream>
#include <chrono>
#include <ctime>
#include <sstream>
#include <iomanip>
#include "time.h"

// Returns current date string shifted by n days ("DD-MM-YYYY")
std::string GetCurrentDateStringFCCode(int n) {
    // Get current time_point from system_clock
    auto now = std::chrono::system_clock::now();
    // Add n days (24 hours per day)
    now += std::chrono::hours(24 * n);

    // Convert to time_t for conversion to tm structure
    std::time_t tt = std::chrono::system_clock::to_time_t(now);

    // Convert to local time (thread-safe alternative: localtime_s)
    std::tm local_tm;
    localtime_s(&local_tm, &tt);  // Corrected order

    // Format the date as "DD-MM-YYYY"
    std::ostringstream oss;
    oss << (local_tm.tm_year + 1900) << "_"
        << std::setfill('0') << std::setw(2) << (local_tm.tm_mon + 1) << "_"
        << std::setfill('0') << std::setw(2) << local_tm.tm_mday ;

    return oss.str();
}

// Returns current date string shifted by n days ("DD-MM-YYYY")
std::string GetCurrentDateString(int n) {
    // Get current time_point from system_clock
    auto now = std::chrono::system_clock::now();
    // Add n days (24 hours per day)
    now += std::chrono::hours(24 * n);

    // Convert to time_t for conversion to tm structure
    std::time_t tt = std::chrono::system_clock::to_time_t(now);

    // Convert to local time (thread-safe alternative: localtime_s)
    std::tm local_tm;
    localtime_s(&local_tm, &tt);  // Corrected order

    // Format the date as "DD-MM-YYYY"
    std::ostringstream oss;
    oss << std::setfill('0') << std::setw(2) << local_tm.tm_mday << "-"
        << std::setfill('0') << std::setw(2) << (local_tm.tm_mon + 1) << "-"
        << (local_tm.tm_year + 1900);

    return oss.str();
}

// Struct to hold GMT time details
struct TimeGMT {
    size_t hours;
    size_t minutes;
    size_t seconds;
    size_t microseconds;
};

// Returns current GMT time with microsecond precision
TimeGMT GetCurrentTimeGMT() {
    TimeGMT result;

    // Get current time with high resolution
    auto now = std::chrono::system_clock::now();
    // Extract seconds portion
    auto now_seconds = std::chrono::time_point_cast<std::chrono::seconds>(now);
    // Calculate microseconds remainder
    auto microsec = std::chrono::duration_cast<std::chrono::microseconds>(now - now_seconds).count();

    // Convert seconds part to time_t
    std::time_t tt = std::chrono::system_clock::to_time_t(now_seconds);

    // Convert to GMT (thread-safe alternative: gmtime_s)
    std::tm gmt_tm;
    gmtime_s(&gmt_tm, &tt);  // Corrected order

    result.hours = static_cast<size_t>(gmt_tm.tm_hour);
    result.minutes = static_cast<size_t>(gmt_tm.tm_min);
    result.seconds = static_cast<size_t>(gmt_tm.tm_sec);
    result.microseconds = static_cast<size_t>(microsec);

    return result;
}

#endif // DATETIME_H
