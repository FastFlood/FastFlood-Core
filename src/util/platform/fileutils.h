#ifndef FILE_UTILITIES_H
#define FILE_UTILITIES_H

#include <string>
#include <vector>
#include <chrono>
#include <functional>
#include <filesystem>
#include <memory>
#include <mutex>
#include <thread>
#include <atomic>
#include <unordered_map>
#include "../../io/ffio.h"
// POSIX/Windows includes:
#if defined(_WIN32)
#  include <windows.h>
#else
#  include <sys/types.h>
#  include <sys/stat.h>
#  include <fcntl.h>
#  include <unistd.h>
#  if defined(__linux__)
#    include <sys/inotify.h>
#  endif
#endif

#ifdef __EMSCRIPTEN__
#  include <emscripten.h>
#endif
inline std::string formatFileSize(uint64_t bytes) {
    const char* units[] = {"B", "KB", "MB", "GB", "TB", "PB", "EB"};
    size_t unitIndex = 0;
    double size = static_cast<double>(bytes);

    while (size >= 1024.0 && unitIndex < sizeof(units) / sizeof(units[0]) - 1) {
        size /= 1024.0;
        ++unitIndex;
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(size < 10 && unitIndex > 0 ? 2 : 0) << size << " " << units[unitIndex];
    return oss.str();
}

    // -------------------------------------------------------------------------
    // Data Structures
    // -------------------------------------------------------------------------

    /// A simple structure holding basic file information.
    struct FileInfo {
        std::string name; ///< File name (without path)
        uintmax_t   size; ///< File size in bytes
        // Last modified time (as a std::chrono time_point)
        std::filesystem::file_time_type lastModified;
    };

    /// A unique identifier for folder watches.
    using FileWatchID = size_t;

    /// Type alias for the callback function called when a folder’s contents change.
    using FolderWatchCallback = std::function<void()>;

    // -------------------------------------------------------------------------
    // GetAllFilesInDir()
    // -------------------------------------------------------------------------

    /**
     * @brief Returns a vector of FileInfo for all regular files in the given folder.
     *
     * Uses std::filesystem to iterate the directory. In case of error, the returned
     * vector might be empty.
     *
     * @param folder The folder path.
     * @return std::vector<FileInfo> A vector containing information about each file.
     */
    inline std::vector<FileInfo> GetAllFilesInDir(const std::string &folder) {
        std::vector<FileInfo> files;
        namespace fs = std::filesystem;
        try {
            for (const auto &entry : fs::directory_iterator(folder)) {
                if (fs::is_regular_file(entry.path())) {
                    FileInfo info;
                    info.name         = entry.path().filename().string();
                    info.size         = fs::file_size(entry.path());
                    info.lastModified = fs::last_write_time(entry.path());
                    files.push_back(info);
                }
            }
        } catch (const fs::filesystem_error &) {
            // Handle error as needed.
        }
        return files;
    }

    /**
     * @brief Returns a vector of FileInfo for all regular files in the given folder.
     *
     * Uses std::filesystem to iterate the directory. In case of error, the returned
     * vector might be empty.
     *
     * @param folder The folder path.
     * @return std::vector<FileInfo> A vector containing information about each file.
     */
    inline std::vector<FileInfo> GetAllFilesInDirMatching(const std::string &folder,const std::string& wildcardPattern) {
        std::string regexPattern = wildcardToRegex(wildcardPattern);
        std::regex re(regexPattern, std::regex::icase); // case-insensitive match

        std::vector<FileInfo> files;
        namespace fs = std::filesystem;
        try {
            for (const auto &entry : fs::directory_iterator(folder)) {
                std::string filename = entry.path().filename().string();
                if (std::regex_match(filename, re)) {

                        if (fs::is_regular_file(entry.path())) {
                        FileInfo info;
                        info.name         = entry.path().filename().string();
                        info.size         = fs::file_size(entry.path());
                        info.lastModified = fs::last_write_time(entry.path());
                        files.push_back(info);
                    }
                }
            }
        } catch (const fs::filesystem_error &) {
            // Handle error as needed.
        }
        return files;
    }


    // -------------------------------------------------------------------------
    // Folder Watch (RegisterFolderWatch / RemoveFolderWatch)
    // -------------------------------------------------------------------------

    /*
       Folder watching is not part of the standard library, so we provide a
       simple cross–platform solution using OS–specific code (or polling as a fallback).
       For Emscripten, these functions are stubbed out.
    */

#if defined(__EMSCRIPTEN__)
    inline FileWatchID RegisterFolderWatch(const std::string & /*folder*/, FolderWatchCallback /*callback*/) {
        return 0;
    }
    inline void RemoveFolderWatch(FileWatchID /*id*/) {}
#elif defined(_WIN32)
    namespace {
        struct WatchInfo {
            std::string folder;
            FolderWatchCallback callback;
            std::atomic<bool> stop;
            std::thread watchThread;
        };
        std::unordered_map<FileWatchID, WatchInfo*> g_watchMap;
        std::mutex g_watchMapMutex;
        std::atomic<FileWatchID> g_nextWatchID{1};

        // A simple polling-based watch thread.
        void WatchFolderThread(WatchInfo* info) {
            auto lastFiles = GetAllFilesInDir(info->folder);
            while (!info->stop) {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                auto currentFiles = GetAllFilesInDir(info->folder);
                if (currentFiles.size() != lastFiles.size()) {
                    info->callback();
                } else {
                    // Compare file names, sizes and modification times.
                    for (size_t i = 0; i < currentFiles.size(); i++) {
                        if (currentFiles[i].name != lastFiles[i].name ||
                            currentFiles[i].size != lastFiles[i].size ||
                            currentFiles[i].lastModified != lastFiles[i].lastModified) {
                            info->callback();
                            break;
                        }
                    }
                }
                lastFiles = currentFiles;
            }
        }
    } // anonymous namespace

    inline FileWatchID RegisterFolderWatch(const std::string &folder, FolderWatchCallback callback) {
        auto* info = new WatchInfo();
        info->folder = folder;
        info->callback = callback;
        info->stop = false;
        info->watchThread = std::thread(WatchFolderThread, info);
        FileWatchID id = g_nextWatchID++;
        {
            std::lock_guard<std::mutex> lock(g_watchMapMutex);
            g_watchMap[id] = info;
        }
        return id;
    }

    inline void RemoveFolderWatch(FileWatchID id) {
        std::lock_guard<std::mutex> lock(g_watchMapMutex);
        auto it = g_watchMap.find(id);
        if (it != g_watchMap.end()) {
            it->second->stop = true;
            if (it->second->watchThread.joinable())
                it->second->watchThread.join();
            delete it->second;
            g_watchMap.erase(it);
        }
    }
#elif defined(__APPLE__) || defined(__linux__)
    namespace {
        struct WatchInfo {
            std::string folder;
            FolderWatchCallback callback;
            std::atomic<bool> stop;
            std::thread watchThread;
#if defined(__linux__)
            int inotifyFd;
            int watchDescriptor;
#endif
        };
        std::unordered_map<FileWatchID, WatchInfo*> g_watchMap;
        std::mutex g_watchMapMutex;
        std::atomic<FileWatchID> g_nextWatchID{1};

#if defined(__linux__)
        void WatchFolderThread(WatchInfo* info) {
            char buffer[1024];
            while (!info->stop) {
                int length = read(info->inotifyFd, buffer, sizeof(buffer));
                if (length > 0) {
                    info->callback();
                }
                std::this_thread::sleep_for(std::chrono::milliseconds(500));
            }
        }
#else
        // macOS fallback: use polling.
        void WatchFolderThread(WatchInfo* info) {
            auto lastFiles = GetAllFilesInDir(info->folder);
            while (!info->stop) {
                std::this_thread::sleep_for(std::chrono::seconds(1));
                auto currentFiles = GetAllFilesInDir(info->folder);
                if (currentFiles.size() != lastFiles.size()) {
                    info->callback();
                } else {
                    for (size_t i = 0; i < currentFiles.size(); i++) {
                        if (currentFiles[i].name != lastFiles[i].name ||
                            currentFiles[i].size != lastFiles[i].size ||
                            currentFiles[i].lastModified != lastFiles[i].lastModified) {
                            info->callback();
                            break;
                        }
                    }
                }
                lastFiles = currentFiles;
            }
        }
#endif
    } // anonymous namespace

    inline FileWatchID RegisterFolderWatch(const std::string &folder, FolderWatchCallback callback) {
        auto* info = new WatchInfo();
        info->folder = folder;
        info->callback = callback;
        info->stop = false;
#if defined(__linux__)
        info->inotifyFd = inotify_init();
        if (info->inotifyFd < 0) {
            delete info;
            return 0;
        }
        info->watchDescriptor = inotify_add_watch(info->inotifyFd, folder.c_str(), IN_MODIFY | IN_CREATE | IN_DELETE);
#endif
        info->watchThread = std::thread(WatchFolderThread, info);
        FileWatchID id = g_nextWatchID++;
        {
            std::lock_guard<std::mutex> lock(g_watchMapMutex);
            g_watchMap[id] = info;
        }
        return id;
    }

    inline void RemoveFolderWatch(FileWatchID id) {
        std::lock_guard<std::mutex> lock(g_watchMapMutex);
        auto it = g_watchMap.find(id);
        if (it != g_watchMap.end()) {
            it->second->stop = true;
            if (it->second->watchThread.joinable())
                it->second->watchThread.join();
#if defined(__linux__)
            if (it->second->watchDescriptor >= 0)
                inotify_rm_watch(it->second->inotifyFd, it->second->watchDescriptor);
            if (it->second->inotifyFd >= 0)
                close(it->second->inotifyFd);
#endif
            delete it->second;
            g_watchMap.erase(it);
        }
    }
#else
    // Fallback stubs if no watch mechanism is available.
    inline FileWatchID RegisterFolderWatch(const std::string &, FolderWatchCallback) {
        return 0;
    }
    inline void RemoveFolderWatch(FileWatchID) {}
#endif

    // -------------------------------------------------------------------------
    // File Locking (LockFile)
    // -------------------------------------------------------------------------

    /**
     * @brief RAII–style file lock.
     *
     * The FileLock object holds an OS–specific handle (or file descriptor) on the file
     * which was opened in read–only mode. When the FileLock object is destroyed, the lock
     * is released automatically.
     */
    struct FileLock {
    public:
#ifdef _WIN32
        HANDLE hFile = INVALID_HANDLE_VALUE;
#else
        int fd = -1;
#endif
        ~FileLock() {
#ifdef _WIN32
            if (hFile != INVALID_HANDLE_VALUE) {
                OVERLAPPED ov = {0};
                UnlockFileEx(hFile, 0, MAXDWORD, MAXDWORD, &ov);
                CloseHandle(hFile);
                hFile = INVALID_HANDLE_VALUE;
            }
#else
            if (fd != -1) {
                // Releasing the lock by closing the file descriptor.
                close(fd);
                fd = -1;
            }
#endif
        }
    };

    /**
     * @brief Locks a file so that it is exclusively in use by the application.
     *
     * This function opens the file in read–only mode (so as not to alter its contents)
     * and acquires an exclusive lock on the entire file. It returns a shared_ptr to a
     * FileLock object that maintains the lock; when the object is destroyed the lock
     * is automatically released.
     *
     * @param filepath The path to the file to lock.
     * @return std::shared_ptr<FileLock> The RAII–style file lock, or nullptr on failure.
     */
    inline std::shared_ptr<FileLock> LockFile(const std::string &filepath) {
#ifdef __EMSCRIPTEN__
        // Emscripten does not support file locking.
        return nullptr;
#elif defined(_WIN32)
        auto lock = std::make_shared<FileLock>();
        // Open the file in read–only mode so as not to alter its contents.
        lock->hFile = CreateFileA(filepath.c_str(), GENERIC_READ, 0, NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
        if (lock->hFile == INVALID_HANDLE_VALUE)
            return nullptr;
        OVERLAPPED ov = {0};
        if (!LockFileEx(lock->hFile, LOCKFILE_EXCLUSIVE_LOCK, 0, MAXDWORD, MAXDWORD, &ov)) {
            CloseHandle(lock->hFile);
            return nullptr;
        }
        return lock;
#else
        auto lock = std::make_shared<FileLock>();
        // Open the file in read–only mode.
        lock->fd = open(filepath.c_str(), O_RDONLY);
        if (lock->fd == -1)
            return nullptr;
        struct flock fl;
        fl.l_type   = F_WRLCK;   // Request an exclusive lock.
        fl.l_whence = SEEK_SET;
        fl.l_start  = 0;
        fl.l_len    = 0;         // Lock the entire file.
        if (fcntl(lock->fd, F_SETLK, &fl) == -1) {
            close(lock->fd);
            return nullptr;
        }
        return lock;
#endif
    }

#endif // FILE_UTILITIES_H
