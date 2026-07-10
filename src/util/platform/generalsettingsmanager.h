#ifndef GENERAL_SETTINGS_MANAGER_H
#define GENERAL_SETTINGS_MANAGER_H

/*
 *  GeneralSettingsManager
 *  ----------------------
 *  Thread-safe, header-only settings database.
 *  - Persists scalar types and std::vector<> types.
 *  - Touches the settings file only when a value actually changed.
 *
 *  Requires a tiny helper:  getExecutableDirectory()  (see executablepath.h)
 */

#include "executablepath.h"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <condition_variable>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <vector>

/* -------------------------------------------------------------------------- */
/*  helpers (stay internal to header)                                         */
/* -------------------------------------------------------------------------- */
namespace gsm_detail
{
    /* ---- lexical_cast-lite ---------------------------------------------- */
    template <class T>
    inline T fromString(const std::string& token)
    {
        std::istringstream iss(token);
        T value{};
        iss >> value;
        return value;
    }
    template <>
    inline std::string fromString<std::string>(const std::string& token)
    {
        return token;
    }

    /* ---- split / join ---------------------------------------------------- */
    template <class T>
    inline std::vector<T> split(const std::string& s, char delim = ',')
    {
        std::vector<T> out;
        std::string token;
        std::istringstream iss(s);
        while (std::getline(iss, token, delim))
            out.emplace_back(fromString<T>(token));
        return out;
    }

    template <class T>
    inline std::string join(const std::vector<T>& v, char delim = ',')
    {
        std::ostringstream oss;
        for (size_t i = 0; i < v.size(); ++i)
        {
            if (i) oss << delim;
            oss << v[i];
        }
        return oss.str();
    }
}  // namespace gsm_detail

/* -------------------------------------------------------------------------- */
/*  GeneralSettingsManager                                                    */
/* -------------------------------------------------------------------------- */
class GeneralSettingsManager
{
public:
    static GeneralSettingsManager& getInstance()
    {
        static GeneralSettingsManager instance;
        return instance;
    }

    template <typename T>
    T getSetting(const std::string& key, const T& defaultValue);

    template <typename T>
    void setSetting(const std::string& key, const T& value);

private:
    /* scalar maps */
    std::unordered_map<std::string, int>               intSettings_;
    std::unordered_map<std::string, bool>              boolSettings_;
    std::unordered_map<std::string, float>             floatSettings_;
    std::unordered_map<std::string, double>            doubleSettings_;
    std::unordered_map<std::string, size_t>            sizeTSettings_;
    std::unordered_map<std::string, std::string>       stringSettings_;

    /* vector maps */
    std::unordered_map<std::string, std::vector<int>>       intVectorSettings_;
    std::unordered_map<std::string, std::vector<float>>     floatVectorSettings_;
    std::unordered_map<std::string, std::vector<size_t>>    sizeTVectorSettings_;
    std::unordered_map<std::string, std::vector<std::string>> stringVectorSettings_;

    /* threading */
    std::mutex              mutex_;
    std::condition_variable cv_;
    bool                    pendingWrites_{false};
    bool                    running_{true};
    std::thread             writerThread_;

    /* misc */
    std::string settingsFile_;

    /* ctor / dtor */
    GeneralSettingsManager();
    ~GeneralSettingsManager();

    /* IO helpers */
    void loadSettings();
    void saveSettings();
    void writerTask();

    /* non-copyable */
    GeneralSettingsManager(const GeneralSettingsManager&)            = delete;
    GeneralSettingsManager& operator=(const GeneralSettingsManager&) = delete;
};

/* -------------------------------------------------------------------------- */
/*  implementation                                                            */
/* -------------------------------------------------------------------------- */
inline GeneralSettingsManager::GeneralSettingsManager()
{
    settingsFile_ = getExecutableDirectory() + "/settings.ini";
    loadSettings();
    writerThread_ = std::thread(&GeneralSettingsManager::writerTask, this);
}

inline GeneralSettingsManager::~GeneralSettingsManager()
{
    {
        std::lock_guard<std::mutex> lock(mutex_);
        running_ = false;
        cv_.notify_one();
    }
    if (writerThread_.joinable()) writerThread_.join();
    saveSettings();   // final flush (if still dirty)
}

/* ---- loadSettings -------------------------------------------------------- */
inline void GeneralSettingsManager::loadSettings()
{
    std::ifstream file(settingsFile_);
    if (!file.is_open()) return;

    std::string key, value;
    while (file >> key >> value)
    {
        if (key.empty()) continue;

        /* vector value?  we encode vectors as "1,2,3" */
        if (value.find(',') != std::string::npos)
        {
            /* decide whether numeric or string vector */
            const auto numericOnly =
                std::all_of(value.begin(), value.end(),
                            [](char c)
                            { return std::isdigit(c) || c == ',' || c == '-' || c == '+' || c == '.'; });

            if (numericOnly)
            {
                if (value.find('.') != std::string::npos)
                {   /* floats */
                    floatVectorSettings_[key]  = gsm_detail::split<float>(value);
                }
                else
                {   /* integers */
                    intVectorSettings_[key]    = gsm_detail::split<int>(value);
                    sizeTVectorSettings_[key]  = gsm_detail::split<size_t>(value);
                }
            }
            else
            {   /* plain strings */
                stringVectorSettings_[key] = gsm_detail::split<std::string>(value);
            }
            continue;
        }

        /* scalar values */
        if (value == "true" || value == "false")
        {
            boolSettings_[key] = (value == "true");
        }
        else if (std::isdigit(value[0]) || value[0] == '-' || value[0] == '+')
        {
            if (value.find('.') != std::string::npos)
                floatSettings_[key] = std::stof(value);
            else
                intSettings_[key]   = std::stoi(value);
        }
        else
        {
            stringSettings_[key] = value;
        }
    }
}

/* ---- saveSettings -------------------------------------------------------- */
inline void GeneralSettingsManager::saveSettings()
{
    std::ofstream file(settingsFile_, std::ios::trunc);
    if (!file.is_open()) return;

    /* scalars */
    for (const auto& [k, v] : intSettings_)    file << k << " " << v << "\n";
    for (const auto& [k, v] : boolSettings_)   file << k << " " << (v ? "true" : "false") << "\n";
    for (const auto& [k, v] : floatSettings_)  file << k << " " << v << "\n";
    for (const auto& [k, v] : doubleSettings_) file << k << " " << v << "\n";
    for (const auto& [k, v] : sizeTSettings_)  file << k << " " << v << "\n";
    for (const auto& [k, v] : stringSettings_) file << k << " " << v << "\n";

    /* vectors */
    for (const auto& [k, v] : intVectorSettings_)
        file << k << " " << gsm_detail::join(v) << "\n";

    for (const auto& [k, v] : floatVectorSettings_)
        file << k << " " << gsm_detail::join(v) << "\n";

    for (const auto& [k, v] : sizeTVectorSettings_)
        file << k << " " << gsm_detail::join(v) << "\n";

    for (const auto& [k, v] : stringVectorSettings_)
        file << k << " " << gsm_detail::join(v) << "\n";
}

/* ---- background writer thread ------------------------------------------- */
inline void GeneralSettingsManager::writerTask()
{
    std::unique_lock<std::mutex> lock(mutex_);
    while (running_)
    {
        cv_.wait(lock, [this] { return pendingWrites_ || !running_; });
        if (!running_) break;

        saveSettings();
        pendingWrites_ = false;
    }
}

/* ---- getSetting ---------------------------------------------------------- */
template <typename T>
inline T GeneralSettingsManager::getSetting(const std::string& key,
                                            const T&           defaultValue)
{
    std::lock_guard<std::mutex> lock(mutex_);

    if constexpr (std::is_same_v<T, int>)
        return intSettings_.count(key) ? intSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, bool>)
        return boolSettings_.count(key) ? boolSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, float>)
        return floatSettings_.count(key) ? floatSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, double>)
        return doubleSettings_.count(key) ? doubleSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, size_t>)
        return sizeTSettings_.count(key) ? sizeTSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, std::string>)
        return stringSettings_.count(key) ? stringSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, std::vector<int>>)
        return intVectorSettings_.count(key) ? intVectorSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, std::vector<float>>)
        return floatVectorSettings_.count(key) ? floatVectorSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, std::vector<size_t>>)
        return sizeTVectorSettings_.count(key) ? sizeTVectorSettings_[key] : defaultValue;
    else if constexpr (std::is_same_v<T, std::vector<std::string>>)
        return stringVectorSettings_.count(key) ? stringVectorSettings_[key] : defaultValue;
    else
        return defaultValue;  // unsupported T
}

/* ---- setSetting ---------------------------------------------------------- */
template <typename T>
inline void GeneralSettingsManager::setSetting(const std::string& key,
                                               const T&           value)
{
    std::lock_guard<std::mutex> lock(mutex_);

    auto flagIfChanged = [this](bool changed)
    {
        if (changed)
        {
            pendingWrites_ = true;
            cv_.notify_one();
        }
    };

    if constexpr (std::is_same_v<T, int>)
    {
        bool changed = (intSettings_[key] != value);
        intSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, bool>)
    {
        bool changed = (boolSettings_[key] != value);
        boolSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, float>)
    {
        bool changed = (floatSettings_[key] != value);
        floatSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, double>)
    {
        bool changed = (doubleSettings_[key] != value);
        doubleSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, size_t>)
    {
        bool changed = (sizeTSettings_[key] != value);
        sizeTSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, std::string>)
    {
        bool changed = (stringSettings_[key] != value);
        stringSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, std::vector<int>>)
    {
        bool changed = (intVectorSettings_[key] != value);
        intVectorSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, std::vector<float>>)
    {
        bool changed = (floatVectorSettings_[key] != value);
        floatVectorSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, std::vector<size_t>>)
    {
        bool changed = (sizeTVectorSettings_[key] != value);
        sizeTVectorSettings_[key] = value;
        flagIfChanged(changed);
    }
    else if constexpr (std::is_same_v<T, std::vector<std::string>>)
    {
        bool changed = (stringVectorSettings_[key] != value);
        stringVectorSettings_[key] = value;
        flagIfChanged(changed);
    }
    /* else: unsupported T – compile-time error is fine */
}

#endif  // GENERAL_SETTINGS_MANAGER_H
