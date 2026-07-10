#ifndef PORTABLE_FILE_DIALOG_H
#define PORTABLE_FILE_DIALOG_H

#include <string>
#include <vector>
#include <stdexcept>
#include <mutex>
#include <condition_variable>
#include <algorithm>
#include <sstream>

//
// This header provides two functions in the PortableDialog namespace:
//
//   • openFileDialog(filetypes, save, multiple)
//   • openFolderDialog(save, multiple)
//
// The “filetypes” parameter is a vector of (extension, description) pairs
// (e.g. { {"png", "PNG Files"}, {"jpg", "JPEG Files"} }).
//
// The file dialog supports both open– and save–modes (save mode shows an
// overwrite prompt and disallows multiple file selection).
//
// The folder dialog now supports a “save” option (default true) to indicate
// that the dialog is used to choose a destination folder (if the native API
// supports such a distinction). Note that on many platforms this flag is ignored,
// and multiple–folder selection is not supported.
//
// In all cases the “last used directory” is stored (thread–safely) so that subsequent
// dialogs can initialize in that folder.
//

/// Open a file–selection dialog.
/// @param filetypes a vector of (extension, description) pairs
/// @param save if true, a “save” dialog is shown (multiple must be false)
/// @param multiple if true, multiple files may be selected (only in open mode)
std::vector<std::string> openFileDialog(
    const std::vector<std::pair<std::string, std::string>>& filetypes,
    bool save = false,
    bool multiple = false);

/// Open a folder–selection dialog.
/// @param save if true (the default), indicates the dialog is used for saving;
///             note that on many platforms this flag is ignored.
/// @param multiple if true, multiple folders may be selected (not supported on most platforms)
std::vector<std::string> openFolderDialog(bool save = true, bool multiple = false);


// ─────────────────────────────────────────────────────────────────────────────
//  IMPLEMENTATION
// ─────────────────────────────────────────────────────────────────────────────

// A global (thread–safe) “last directory” is stored so that dialogs can start
// in the last–used location.
static std::mutex lastDirMutex;
static std::string lastDirectory; // e.g., "C:\\Users\\User\\Documents" or "/home/user"

static void updateLastDirectory(const std::string &path) {
    std::lock_guard<std::mutex> lock(lastDirMutex);
    // Extract the directory portion (works for both Windows and Unix paths)
    size_t pos = path.find_last_of("\\/");
    if (pos != std::string::npos)
        lastDirectory = path.substr(0, pos);
}

//
// ─── WINDOWS IMPLEMENTATION ─────────────────────────────────────────────────────
//
#ifdef _WIN32

#include <windows.h>
#include <commdlg.h>
#include <shlobj.h> // for folder selection

inline std::vector<std::string> openFileDialog(
    const std::vector<std::pair<std::string, std::string>>& filetypes,
    bool save,
    bool multiple)
{
    constexpr int FILE_DIALOG_MAX_BUFFER = 16384;
    if (save && multiple)
        throw std::invalid_argument("Cannot have both save and multiple=true.");

    std::vector<std::string> result;
    OPENFILENAMEW ofn;
    ZeroMemory(&ofn, sizeof(ofn));
    wchar_t fileBuffer[FILE_DIALOG_MAX_BUFFER] = {0};
    ofn.lStructSize = sizeof(OPENFILENAMEW);
    ofn.lpstrFile = fileBuffer;
    ofn.nMaxFile = FILE_DIALOG_MAX_BUFFER;
    ofn.nFilterIndex = 1;

    // Build the filter string.
    std::string filter;
    if (!save && filetypes.size() > 1) {
        filter.append("Supported Files (");
        for (size_t i = 0; i < filetypes.size(); ++i) {
            filter.append("*.");
            filter.append(filetypes[i].first);
            if (i + 1 < filetypes.size())
                filter.append(";");
        }
        filter.append(")");
        filter.push_back('\0');
        for (size_t i = 0; i < filetypes.size(); ++i) {
            filter.append("*.");
            filter.append(filetypes[i].first);
            if (i + 1 < filetypes.size())
                filter.append(";");
        }
        filter.push_back('\0');
    }
    for (auto &pair : filetypes) {
        filter.append(pair.second);
        filter.append(" (*.");
        filter.append(pair.first);
        filter.append(")");
        filter.push_back('\0');
        filter.append("*.");
        filter.append(pair.first);
        filter.push_back('\0');
    }
    filter.push_back('\0');

    // Convert filter to a wide string.
    int sizeNeeded = MultiByteToWideChar(CP_UTF8, 0, filter.c_str(), -1, NULL, 0);
    std::wstring wfilter(sizeNeeded, L'\0');
    MultiByteToWideChar(CP_UTF8, 0, filter.c_str(), -1, &wfilter[0], sizeNeeded);
    ofn.lpstrFilter = wfilter.c_str();

    // Set the initial directory from the global variable, if available.
    {
        std::lock_guard<std::mutex> lock(lastDirMutex);
        if (!lastDirectory.empty()) {
            std::wstring wdir(lastDirectory.begin(), lastDirectory.end());
            ofn.lpstrInitialDir = wdir.c_str();
        }
    }

    if (save) {
        ofn.Flags = OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_OVERWRITEPROMPT;
        if (GetSaveFileNameW(&ofn) == FALSE)
            return {};
    } else {
        ofn.Flags = OFN_EXPLORER | OFN_PATHMUSTEXIST | OFN_FILEMUSTEXIST;
        if (multiple)
            ofn.Flags |= OFN_ALLOWMULTISELECT;
        if (GetOpenFileNameW(&ofn) == FALSE)
            return {};
    }

    // Parse the returned buffer.
    size_t i = 0;
    while (fileBuffer[i] != L'\0') {
        std::wstring wstr(&fileBuffer[i]);
        int utf8Size = WideCharToMultiByte(CP_UTF8, 0, wstr.c_str(), -1, NULL, 0, NULL, NULL);
        std::string filename(utf8Size, '\0');
        WideCharToMultiByte(CP_UTF8, 0, wstr.c_str(), -1, &filename[0], utf8Size, NULL, NULL);
        if (!filename.empty() && filename.back() == '\0')
            filename.pop_back();
        result.push_back(filename);
        i += wstr.size() + 1;
    }
    // When multiple files are selected, the first element is the folder.
    if (result.size() > 1) {
        std::string folder = result[0];
        for (size_t j = 1; j < result.size(); ++j)
            result[j] = folder + "\\" + result[j];
        result.erase(result.begin());
    }
    if (!result.empty())
        updateLastDirectory(result.front());
    return result;
}

inline std::vector<std::string> openFolderDialog(bool /*save*/, bool multiple)
{
    if (multiple)
        throw std::invalid_argument("Multiple folder selection is not supported on Windows.");
    // For folder selection on Windows, we use SHBrowseForFolderW.
    // The 'save' flag is ignored here.
    BROWSEINFOW bi = {0};
    bi.ulFlags = BIF_RETURNONLYFSDIRS | BIF_NEWDIALOGSTYLE;
    wchar_t path[MAX_PATH] = {0};
    bi.lpszTitle = L"Select Folder";
    bi.pszDisplayName = path;
    LPITEMIDLIST pidl = SHBrowseForFolderW(&bi);
    if (pidl != nullptr) {
        if (SHGetPathFromIDListW(pidl, path)) {
            CoTaskMemFree(pidl);
            int utf8Size = WideCharToMultiByte(CP_UTF8, 0, path, -1, NULL, 0, NULL, NULL);
            std::string folder(utf8Size, '\0');
            WideCharToMultiByte(CP_UTF8, 0, path, -1, &folder[0], utf8Size, NULL, NULL);
            if (!folder.empty() && folder.back() == '\0')
                folder.pop_back();
            updateLastDirectory(folder);
            return { folder };
        }
        CoTaskMemFree(pidl);
    }
    return {};
}

//
// ─── macOS IMPLEMENTATION ────────────────────────────────────────────────────────
//
#elif defined(__APPLE__)

#include <cstdio>
#include <array>

std::vector<std::string> openFileDialog(
    const std::vector<std::pair<std::string, std::string>>& filetypes,
    bool save,
    bool multiple)
{
    std::string script;
    if (save) {
        // "Save As" file dialog.
        script = "tell application \"System Events\" to return POSIX path of (choose file name)";
    } else {
        script = "tell application \"System Events\" to return POSIX path of (choose file";
        if (multiple)
            script += " with multiple selections allowed";
        if (!filetypes.empty()) {
            script += " of type {";
            for (size_t i = 0; i < filetypes.size(); ++i) {
                script += "\"" + filetypes[i].first + "\"";
                if (i + 1 < filetypes.size())
                    script += ", ";
            }
            script += "}";
        }
        script += ")";
    }
    std::string command = "osascript -e '" + script + "'";
    std::array<char, 4096> buffer;
    std::string output;
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe)
        return {};
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr)
        output += buffer.data();
    pclose(pipe);
    output.erase(std::remove(output.begin(), output.end(), '\n'), output.end());
    std::vector<std::string> result;
    if (!output.empty()) {
        result.push_back(output);
        updateLastDirectory(output);
    }
    return result;
}

std::vector<std::string> openFolderDialog(bool /*save*/, bool multiple)
{
    if (multiple)
        throw std::invalid_argument("Multiple folder selection is not supported on macOS.");
    // AppleScript does not differentiate between save and open for folder selection.
    std::string script = "tell application \"System Events\" to return POSIX path of (choose folder)";
    std::string command = "osascript -e '" + script + "'";
    std::array<char, 4096> buffer;
    std::string output;
    FILE* pipe = popen(command.c_str(), "r");
    if (!pipe)
        return {};
    while (fgets(buffer.data(), buffer.size(), pipe) != nullptr)
        output += buffer.data();
    pclose(pipe);
    output.erase(std::remove(output.begin(), output.end(), '\n'), output.end());
    std::vector<std::string> result;
    if (!output.empty()) {
        result.push_back(output);
        updateLastDirectory(output);
    }
    return result;
}

//
// ─── EMSCRIPTEN (WebAssembly) IMPLEMENTATION ─────────────────────────────────────
//
#elif defined(__EMSCRIPTEN__)

#include <emscripten.h>
#include <emscripten/html5.h>
#include <thread>
#include <chrono>
#include <atomic>

namespace {
    std::mutex emsMutex;
    std::condition_variable emsCV;
    std::vector<std::string> emsResult;
    bool emsDone = false;
}

extern "C" {
    // Called from JavaScript when file/folder selection is complete.
    void fileDialogCallback(const char* result) {
        std::unique_lock<std::mutex> lock(emsMutex);
        emsResult.clear();
        if (result && result[0] != '\0') {
            std::istringstream iss(result);
            std::string token;
            while (std::getline(iss, token, '\n')) {
                if (!token.empty())
                    emsResult.push_back(token);
            }
        }
        emsDone = true;
        emsCV.notify_one();
    }
}

// This helper creates and clicks an invisible file input element. It accepts
// parameters to enable folder selection (via webkitdirectory), multiple selection,
// and a file filter (a comma–separated list).
// (Browsers do not support a “save–file” dialog so the save flag is ignored.)
static std::vector<std::string> runEmscriptenFileDialog(bool folder, bool multiple, const std::string &accept) {
    {
        std::unique_lock<std::mutex> lock(emsMutex);
        emsDone = false;
    }
    EM_ASM_({
        var multiple = $0;
        var folder = $1;
        var accept = UTF8ToString($2);
        var input = document.createElement('input');
        input.type = 'file';
        if (multiple)
            input.multiple = true;
        if (folder)
            input.webkitdirectory = true;
        if (accept && accept.length > 0)
            input.accept = accept;
        input.style.display = 'none';
        document.body.appendChild(input);
        input.onchange = function(e) {
            var files = input.files;
            var result = "";
            for (var i = 0; i < files.length; i++) {
                result += files[i].name;
                if (i < files.length - 1)
                    result += "\n";
            }
            Module.ccall('fileDialogCallback', 'void', ['string'], [result]);
            document.body.removeChild(input);
        };
        input.click();
    }, multiple ? 1 : 0, folder ? 1 : 0, accept.c_str());

    std::unique_lock<std::mutex> lock(emsMutex);
    // Wait for the JS callback. (This requires compiling with Asyncify enabled.)
    while (!emsDone)
        emscripten_sleep(10);
    return emsResult;
}

std::vector<std::string> openFileDialog(
    const std::vector<std::pair<std::string, std::string>>& filetypes,
    bool save,
    bool multiple)
{
    (void)save; // "save" is ignored in browsers.
    std::string accept;
    for (const auto &p : filetypes) {
        if (!accept.empty())
            accept += ",";
        accept += "." + p.first;
    }
    return runEmscriptenFileDialog(false, multiple, accept);
}

std::vector<std::string> openFolderDialog(bool /*save*/, bool multiple)
{
    if (multiple)
        throw std::invalid_argument("Multiple folder selection is not supported in browsers.");
    // The save flag is ignored for folder dialogs in browsers.
    return runEmscriptenFileDialog(true, false, "");
}

//
// ─── FALLBACK (Linux / Other using Zenity) ──────────────────────────────────────
//
#else

#include <cstdio>
#include <array>

std::vector<std::string> openFileDialog(
    const std::vector<std::pair<std::string, std::string>>& filetypes,
    bool save,
    bool multiple)
{
    constexpr int BUFFER_SIZE = 16384;
    char buffer[BUFFER_SIZE] = {0};
    std::string cmd = "zenity --file-selection ";
    if (multiple)
        cmd += "--multiple --separator=\"\n\" ";
    if (save)
        cmd += "--save ";
    if (!filetypes.empty()) {
        cmd += "--file-filter=\"";
        for (size_t i = 0; i < filetypes.size(); ++i) {
            cmd += "*." + filetypes[i].first;
            if (i + 1 < filetypes.size())
                cmd += " ";
        }
        cmd += "\" ";
    }
    {
        std::lock_guard<std::mutex> lock(lastDirMutex);
        if (!lastDirectory.empty())
            cmd += "--filename=\"" + lastDirectory + "/\" ";
    }
    FILE *pipe = popen(cmd.c_str(), "r");
    if (!pipe)
        throw std::runtime_error("popen() failed for zenity!");
    std::string output;
    while (fgets(buffer, BUFFER_SIZE, pipe) != nullptr)
        output += buffer;
    pclose(pipe);
    if (!output.empty() && output.back() == '\n')
        output.pop_back();

    std::vector<std::string> result;
    if (multiple) {
        size_t start = 0;
        size_t pos = 0;
        while ((pos = output.find('\n', start)) != std::string::npos) {
            result.push_back(output.substr(start, pos - start));
            start = pos + 1;
        }
        if (start < output.size())
            result.push_back(output.substr(start));
    } else {
        result.push_back(output);
    }
    if (!result.empty())
        updateLastDirectory(result.front());
    return result;
}

std::vector<std::string> openFolderDialog(bool save, bool multiple)
{
    if (multiple)
        throw std::invalid_argument("Multiple folder selection is not supported with zenity.");
    constexpr int BUFFER_SIZE = 16384;
    char buffer[BUFFER_SIZE] = {0};
    std::string cmd = "zenity --file-selection --directory ";
    if (save)
        cmd += "--save ";
    {
        std::lock_guard<std::mutex> lock(lastDirMutex);
        if (!lastDirectory.empty())
            cmd += "--filename=\"" + lastDirectory + "/\" ";
    }
    FILE *pipe = popen(cmd.c_str(), "r");
    if (!pipe)
        throw std::runtime_error("popen() failed for zenity!");
    std::string output;
    while (fgets(buffer, BUFFER_SIZE, pipe) != nullptr)
        output += buffer;
    pclose(pipe);
    if (!output.empty() && output.back() == '\n')
        output.pop_back();
    std::vector<std::string> result;
    if (!output.empty()) {
        result.push_back(output);
        updateLastDirectory(output);
    }
    return result;
}

#endif // OS-specific implementations


#endif // PORTABLE_FILE_DIALOG_H
