#ifndef EXECUTABLEPATH_H
#define EXECUTABLEPATH_H
#include <iostream>
#include <string>

#if defined(_WIN32)
    #include <windows.h>
    #include <libloaderapi.h>
#elif defined(__APPLE__)
    #include <mach-o/dyld.h>
#elif defined(__linux__)
    #include <unistd.h>
#elif defined(__EMSCRIPTEN__)
    #include <emscripten.h>
#endif

#include <cstdlib>
#include <string>

inline std::string GetHomeDirectory() {
#ifdef _WIN32
    // Try the USERPROFILE environment variable first.
    if (const char* homeDir = std::getenv("USERPROFILE")) {
        return std::string(homeDir);
    }
    // Fallback: combine HOMEDRIVE and HOMEPATH if available.
    const char* drive = std::getenv("HOMEDRIVE");
    const char* path  = std::getenv("HOMEPATH");
    if (drive && path) {
        return std::string(drive) + std::string(path);
    }
    return ""; // Unable to determine home directory.
#else
    // For Unix-like systems (macOS, Linux, etc.)
    if (const char* homeDir = std::getenv("HOME")) {
        return std::string(homeDir);
    }
#ifdef __EMSCRIPTEN__
    // Emscripten may not have the HOME environment variable set.
    // Provide a reasonable default for the simulated filesystem.
    return std::string("/home/web_user");
#else
    return ""; // Unable to determine home directory.
#endif
#endif
}


inline std::string getExecutableDirectory() {
    char path[1024] = {0};

#if defined(_WIN32)
    wchar_t wpath[MAX_PATH];
    if (GetModuleFileNameW(nullptr, wpath, MAX_PATH)) {
        std::wstring wstr(wpath);
        std::string fullPath(wstr.begin(), wstr.end());
        size_t pos = fullPath.find_last_of("\\/");
        return (pos == std::string::npos) ? fullPath : fullPath.substr(0, pos);
    }
#elif defined(__APPLE__)
    uint32_t size = sizeof(path);
    if (_NSGetExecutablePath(path, &size) == 0) {
        std::string fullPath(path);
        size_t pos = fullPath.find_last_of("/");
        return (pos == std::string::npos) ? fullPath : fullPath.substr(0, pos);
    }
#elif defined(__linux__)
    ssize_t len = readlink("/proc/self/exe", path, sizeof(path) - 1);
    if (len != -1) {
        path[len] = '\0';
        std::string fullPath(path);
        size_t pos = fullPath.find_last_of("/");
        return (pos == std::string::npos) ? fullPath : fullPath.substr(0, pos);
    }
#elif defined(__EMSCRIPTEN__)
    return ".";  // WebAssembly runs in a virtual filesystem; no real executable path exists.
#endif

    return "";
}
#endif
