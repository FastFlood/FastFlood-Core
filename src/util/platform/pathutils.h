#ifndef PATH_UTILS_H
#define PATH_UTILS_H

#include <string>
#include <algorithm>

#ifdef _WIN32
#include <cctype>
#endif
namespace PathUtils {


inline std::string GetPathSeperator()
{
#ifdef _WIN32
    return "\\";
#else
    return "/";
#endif
}

inline bool isRelativePath(const std::string& pathStr) {
#ifdef _WIN32
    // On Windows:
    // - Absolute paths often start with a drive letter (e.g., "C:\")
    // - UNC paths start with two slashes (e.g., "\\server\share" or "//server/share")
    if (pathStr.size() >= 2 && std::isalpha(pathStr[0]) && pathStr[1] == ':')
        return false;
    if (pathStr.size() >= 2 &&
        ((pathStr[0] == '\\' && pathStr[1] == '\\') || (pathStr[0] == '/' && pathStr[1] == '/')))
        return false;
#else
    // On Unix-like systems, an absolute path starts with '/'
    if (!pathStr.empty() && pathStr[0] == '/')
        return false;
#endif
    return true;
}

    // Preferred path separator for the current platform.
#ifdef _WIN32
    constexpr char preferred_separator = '\\';
    constexpr char alt_separator = '/';
#else
    constexpr char preferred_separator = '/';
    constexpr char alt_separator = '\\';
#endif

    /// Converts all path separators in the given path to the preferred separator
    /// for the current platform.
    inline std::string NormalizePath(std::string path) {
        std::replace(path.begin(), path.end(), alt_separator, preferred_separator);
        return path;
    }

    /// Returns the drive letter (e.g. "C:") on Windows.
    /// On non-Windows platforms, or if no drive is present, returns an empty string.
    inline std::string GetDrive(const std::string& path) {
#ifdef _WIN32
        // Check for drive letter at the beginning (e.g., "C:")
        if (path.size() >= 2 && std::isalpha(static_cast<unsigned char>(path[0])) && path[1] == ':')
            return path.substr(0, 2);
#endif
        return "";
    }

    /*---------------------------------------------------------------
     *  Return true if the last path component contains an extension.
     *  Extension = a dot ('.') that is not the first character of
     *  the component and is followed by at least one other byte.
     *  The path is UTF-8; scanning bytes is safe because '.' (0x2E)
     *  cannot appear inside a UTF-8 multibyte sequence.
     *--------------------------------------------------------------*/
    inline bool HasExtension(const std::string &path)
    {
        if (path.empty())
            return false;

        /* ignore trailing separators ---------------------------------------*/
        std::size_t end = path.find_last_not_of("/\\");
        if (end == std::string::npos)
        {
            return false;                       /* path is only separators */
        }

        /* locate start of last component -----------------------------------*/
        std::size_t start = path.find_last_of("/\\", end);
        std::size_t comp_begin = (start == std::string::npos) ? 0 : start + 1;

        /* find last '.' before `end` ---------------------------------------*/
        std::size_t dot = path.find_last_of('.', end);

        /* no dot, or dot is before component, or dot is first char, or dot
           is final char  no extension                                    */
        if (dot == std::string::npos           ||   // no dot
            dot < comp_begin + 1               ||   // ".hidden" files
            dot >= end)                             // "filename."
        {
            return false;
        }

        return true;
    }

    /// Returns the directory part of the given path (without the trailing separator).
    /// If no directory separator is found, returns an empty string.
    inline std::string GetDirectory(const std::string& path) {
        auto pos = path.find_last_of("/\\");
        if (pos == std::string::npos)
            return "";
        return path.substr(0, pos);
    }

    /// Returns the file name (with extension) of the given path.
    inline std::string GetFileName(const std::string& path) {
        auto pos = path.find_last_of("/\\");
        if (pos == std::string::npos)
            return path;
        return path.substr(pos + 1);
    }

    /// Returns the extension (including the dot, e.g. ".txt") of the file.
    /// If there is no extension (or the file name starts with a dot), returns an empty string.
    inline std::string GetExtension(const std::string& path) {
        std::string filename = GetFileName(path);
        auto pos = filename.find_last_of('.');
        if (pos == std::string::npos || pos == 0)
            return "";
        return filename.substr(pos);
    }

    /// Returns the bare file name (without extension) of the given path.
    inline std::string GetBaseName(const std::string& path) {
        std::string filename = GetFileName(path);
        auto pos = filename.find_last_of('.');
        if (pos == std::string::npos || pos == 0)
            return filename;
        return filename.substr(0, pos);
    }

    /// Returns a new path with the file's extension replaced by new_ext.
    /// If new_ext is non-empty and does not start with a dot, one is added.
    /// If new_ext is empty, the extension is removed.
    inline std::string SetExtension(const std::string& path, const std::string& new_ext) {
        std::string drive = GetDrive(path);
        std::string dir = GetDirectory(path);
        std::string base = GetBaseName(path);

        std::string ext;
        if (!new_ext.empty()) {
            ext = (new_ext.front() == '.') ? new_ext : ('.' + new_ext);
        }

        std::string newFile = base + ext;
        if (!dir.empty())
            return drive + dir + preferred_separator + newFile;
        else
            return drive + newFile;
    }

    /// Returns a new path with the file name (and its extension) replaced by new_filename.
    inline std::string SetFileName(const std::string& path, const std::string& new_filename) {
        std::string drive = GetDrive(path);
        std::string dir = GetDirectory(path);
        if (!dir.empty())
            return drive + dir + preferred_separator + new_filename;
        else
            return drive + new_filename;
    }

    /// Joins two paths, ensuring that the proper separator is used.
    /// If the second path is absolute (or has a drive letter on Windows),
    /// it is returned after normalization.
    inline std::string JoinPaths(const std::string& path1, const std::string& path2) {
#ifdef _WIN32
        // On Windows, check if path2 has a drive letter or starts with a separator.
        if (path2.size() >= 2 && std::isalpha(static_cast<unsigned char>(path2[0])) && path2[1] == ':')
            return NormalizePath(path2);
        if (!path2.empty() && (path2.front() == '\\' || path2.front() == '/'))
            return NormalizePath(path2);
#else
        // On POSIX systems, if path2 starts with '/', it's absolute.
        if (!path2.empty() && path2.front() == '/')
            return NormalizePath(path2);
#endif
        std::string norm1 = NormalizePath(path1);
        std::string norm2 = NormalizePath(path2);
        if (norm1.empty())
            return norm2;
        if (norm1.back() != preferred_separator)
            norm1.push_back(preferred_separator);
        return norm1 + norm2;
    }

} // namespace PathUtils

#endif // PATH_UTILS_H
