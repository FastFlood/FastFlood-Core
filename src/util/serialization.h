// ============================================================================
//  serialization.hpp ‒ header-only reflective helpers for nlohmann::json
//  2025-05-04
//
//  •  No Meta<> class; uses explicit specialisation via macros.
//  •  Works transparently with STL containers such as
//     std::vector<MyClass>, std::map<std::string,MyClass>, …
//  •  Pointer / non-pointer overloads of to_json() preserved.
//  •  Only *registered* types trigger reflection; everything else
//     falls back to nlohmann-json’s own converters.
// ============================================================================

#ifndef SERIALIZATION_HPP
#define SERIALIZATION_HPP

#include <string>
#include <vector>
#include <map>
#include <mutex>
#include <type_traits>
#include <functional>
#include <nlohmann/json.hpp>
#include "iostream"

class SerializableBase
{
public:
    inline virtual nlohmann::json GetTemplate() const{
        return nlohmann::json();
    }
    inline virtual std::string GetTemplateString() const {
        return  GetTemplate().dump();
        }
    inline virtual nlohmann::json GetJson() const{
        return nlohmann::json();
    }
    inline virtual std::string GetJsonString() const
    {
        return "{}";
    }
    inline virtual bool SetFromJson(const nlohmann::json & j)
    {
        return true;
    }
    inline virtual bool SetFromJsonString(const std::string & jstr)
    {
        return true;
    }
    inline virtual std::string GetTypeNameString() const
    {
        return "SerializableBase";
    }
};

// ────────────────────────────────────────────────────────────────
//  Field metadata
// ────────────────────────────────────────────────────────────────
struct FieldDesc
{
    const char*   name{};
    const char*   type{};
    const char*   description{};
    const char*   defaultRaw{};
    bool          required{false};

    std::function<SerializableBase*(void)> creator;

    std::function<nlohmann::json(const void*)>     getter;
    std::function<bool(void*, const std::string&)> setter;
};

using FieldList = std::vector<FieldDesc>;

// ────────────────────────────────────────────────────────────────
//  Registration infrastructure
// ────────────────────────────────────────────────────────────────
template<typename T> struct is_registered : std::false_type {};

// “forward declarations”—the unspecialised versions must *never* be
// instantiated; if they are, you forgot to register the type.
template<typename T> std::string  GetTypeName()
{
    static_assert(is_registered<T>::value,
                  "Type T has not been registered with "
                  "SERIALIZABLE_FIELDS_BEGIN / END.");
    return {};
}

template<typename T> FieldList    RegisterFields()
{
    static_assert(is_registered<T>::value,
                  "Type T has not been registered with "
                  "SERIALIZABLE_FIELDS_BEGIN / END.");
    return {};
}

// ────────────────────────────────────────────────────────────────
//  Thread-safe, per-type cache of FieldList
// ────────────────────────────────────────────────────────────────
template<typename T>
inline FieldList& fields()
{
    static std::map<std::string, FieldList> cache;
    static std::mutex                       mtx;

    const std::string typeName = GetTypeName<T>();

    std::lock_guard<std::mutex> lock(mtx);
    auto [it, inserted] = cache.try_emplace(typeName);
    if (inserted) it->second = RegisterFields<T>();
    return it->second;
}

// ────────────────────────────────────────────────────────────────
//  Reflection-based (de)serialisation for *registered* types
// ────────────────────────────────────────────────────────────────
template<typename T>
inline nlohmann::json ToJson(const T& obj)
{
    nlohmann::json j = nlohmann::json::object();
    for (auto& f : fields<T>())
        j[f.name] = f.getter(&obj);
    return j;
}

template<typename T>
inline bool FromJson(const nlohmann::json& j, T& obj)
{
    for (auto& f : fields<T>())
    {
        auto it = j.find(f.name);
        if (it != j.end() && !f.setter(&obj, it->dump()))
            return false;
    }
    return true;
}

// ────────────────────────────────────────────────────────────────
//  ADL hooks for nlohmann::json – enabled **only** for registered
//  classes, so STL containers of those classes keep working.
// ────────────────────────────────────────────────────────────────
template<typename T,
         std::enable_if_t<is_registered<T>::value, int> = 0>
inline void to_json(nlohmann::json& j, const T& obj)
{
    j = ToJson(obj);
}

template<typename T,
         std::enable_if_t<is_registered<T>::value, int> = 0>
inline void from_json(const nlohmann::json& j, T& obj)
{
    if (!FromJson(j, obj))
        throw std::runtime_error("from_json failed");
}

// ────────────────────────────────────────────────────────────────
//  Convenience helpers that **return** a json value
// ────────────────────────────────────────────────────────────────

// 1) registered, non-pointer
template<typename T,
         std::enable_if_t< is_registered<T>::value &&
                          !std::is_pointer_v<T>, int> = 0>
inline nlohmann::json to_json(const T& obj)
{
    return ToJson(obj);
}

// 2) NOT registered, non-pointer → fall back to nlohmann::json
template<typename T,
         std::enable_if_t<!is_registered<T>::value &&
                          !std::is_pointer_v<T>, int> = 0>
inline nlohmann::json to_json(const T& obj)
{
    return nlohmann::json(obj);   // relies on nlohmann’s own ADL
}

// 3) pointer to registered
template<typename T,
         std::enable_if_t< std::is_pointer_v<T> &&
                           is_registered<std::remove_pointer_t<T>>::value,
                           int> = 0>
inline nlohmann::json to_json(const T ptr)
{
    if (!ptr) return nullptr;
    return ToJson(*ptr);
}


// ────────────────────────────────────────────────────────────────
//  Optional empty “tag” base class
// ────────────────────────────────────────────────────────────────
template<typename T>
class Serializable : public SerializableBase{ /* empty tag – for user convenience */
public:
    inline nlohmann::json GetTemplate() const override{
        std::cout << "do template " << GetTypeNameString() << std::endl;

        nlohmann::json t = nlohmann::json::object();

        for (auto& f : fields<T>()) {
            //get the sub-template by making an object for it
            if(f.creator)
            {
                SerializableBase * subobj = f.creator();
                nlohmann::json templ = subobj->GetTemplate();
                std::cout << subobj->GetTypeNameString() << std::endl;
                delete subobj;
                t[f.name] = {
                    {"type",        f.type},
                    {"description", f.description},
                    {"required",    f.required},
                    {"default",     f.defaultRaw},
                    {"template",    templ}
                };
            }else
            {
                t[f.name] = {
                    {"type",        f.type},
                    {"description", f.description},
                    {"required",    f.required},
                    {"default",     f.defaultRaw}
                };
            }

        }
        return t;
    }

    std::string GetTemplateString() const override{
        return  GetTemplate().dump();
        }

    inline nlohmann::json GetJson() const override
    {
        return ToJson<T>(*(static_cast<const T*>(this)));
    }
    inline std::string GetJsonString() const override
    {
        return ToJson<T>(*(static_cast<const T*>(this))).dump();
    }
    inline bool SetFromJson(const nlohmann::json & j) override
    {
        FromJson(j,*(static_cast<T*>(this)));
        return true;
    }
    inline bool SetFromJsonString(const std::string & jstr) override
    {
        FromJson(nlohmann::json::parse(jstr),*(static_cast<T*>(this)));
        return true;
    }
    inline std::string GetTypeNameString() const
    {
        return GetTypeName<T>();
    }

};

// ────────────────────────────────────────────────────────────────
//  Registration macros
// ────────────────────────────────────────────────────────────────
#define SERIALIZABLE_FIELDS_BEGIN(CLS)                                      \
    template<> struct is_registered<CLS> : std::true_type {};               \
    template<> inline std::string GetTypeName<CLS>() { return #CLS; }       \
    template<> inline FieldList   RegisterFields<CLS>()                     \
    {                                                                       \
        FieldList _f;

#define SERIALIZABLE_FIELD(CLS, NAME, MEMBER, TYPER, TYPE, DEFRAW, REQ, DESC)\
    {                                                                       \
        FieldDesc f{};                                                      \
        f.name        = NAME;                                               \
        f.type        = TYPE;                                               \
        f.description = DESC;                                               \
        f.defaultRaw  = DEFRAW;                                             \
        f.required    = REQ;                                                \
                                                                            \
        f.getter = [](const void* self)->nlohmann::json                     \
        {                                                                   \
            const auto* obj = static_cast<const CLS*>(self);                \
            return nlohmann::json(obj->MEMBER);                             \
        };                                                                  \
                                                                            \
        f.setter = [](void* self, const std::string& jstr)->bool            \
        {                                                                   \
            auto* obj = static_cast<CLS*>(self);                            \
            try                                                             \
            {                                                               \
                obj->MEMBER = nlohmann::json::parse(jstr).get<TYPER>();     \
                return true;                                                \
            }                                                               \
            catch (...)                                                     \
            { return false; }                                               \
        };                                                                  \
        _f.push_back(std::move(f));                                         \
    }


#define SERIALIZABLE_FIELD_COMPOUND(CLS, NAME, MEMBER, TYPER, TYPE, DEFRAW, REQ, DESC)\
    {                                                                       \
        FieldDesc f{};                                                      \
        f.name        = NAME;                                               \
        f.type        = TYPE;                                               \
        f.description = DESC;                                               \
        f.defaultRaw  = DEFRAW;                                             \
        f.required    = REQ;                                                \
                                                                            \
        f.getter = [](const void* self)->nlohmann::json                     \
        {                                                                   \
            const auto* obj = static_cast<const CLS*>(self);                \
            return nlohmann::json(obj->MEMBER);                             \
        };                                                                  \
                                                                            \
        f.setter = [](void* self, const std::string& jstr)->bool            \
        {                                                                   \
            auto* obj = static_cast<CLS*>(self);                            \
            try                                                             \
            {                                                               \
                obj->MEMBER = nlohmann::json::parse(jstr).get<TYPER>();     \
                return true;                                                \
            }                                                               \
            catch (...)                                                     \
            { return false; }                                               \
        };                                                                  \
        f.creator = [](void)->SerializableBase* { return new TYPER(); } ;      \
        _f.push_back(std::move(f));                                         \
    }

#define SERIALIZABLE_FIELDS_END()   return _f; }

class DynamicSerializable : public SerializableBase
{
public:
    using Json = nlohmann::json;

    struct FieldInfo {
        std::string           name;
        std::string           typeName;
        std::string           description;
        std::string           defaultRaw;
        bool                  required{};
        Json                  nestedTemplate;  // if provided

        void*                 storage    = nullptr;
        std::function<void(void*)>                destructor;
        std::function<Json(void*)>                toJson;
        std::function<bool(void*,const Json&)>    fromJson;
    };

    DynamicSerializable() = default;
    ~DynamicSerializable() { clear(); }

    //───────────────────────────────────────────────────────────────────────────
    //  1) Register any Serializable-derived type by name:
    //     – T must inherit SerializableBase (or Serializable<T>).
    //     – Its GetTypeName<T>() must return the same string.
    //───────────────────────────────────────────────────────────────────────────
    template<typename T>
    static void RegisterSerializableType() {
        const std::string typ = GetTypeName<T>();
        std::lock_guard<std::mutex> lock(_reg_mtx());
        _factories()[typ] = [](){
            // allocate a T() as void*
            return static_cast<void*>(new T());
        };
        _destructors()[typ] = [](void* p){
            delete static_cast<T*>(p);
        };
        _toJsonF()[typ] = [](void* p){
            return static_cast<T*>(p)->GetJson();
        };
        _fromJsonF()[typ] = [](void* p, const Json& j){
            return static_cast<T*>(p)->SetFromJson(j);
        };
    }

    //───────────────────────────────────────────────────────────────────────────
    //  2) Build fields from a JSON‐template string:
    //      • built-in primitives are auto-handled
    //      • anything else must have been RegisterSerializableType<T>()
    //───────────────────────────────────────────────────────────────────────────
    bool SetTemplateString(const std::string& tplStr) {
        clear();
        Json tpl = Json::parse(tplStr);

        for (auto& [key,info] : tpl.items()) {
            std::string typ    = info.value("type","");
            std::string desc   = info.value("description","");
            bool        req    = info.value("required",false);
            std::string defRaw = info.value("default","");
            Json        subTpl = info.value("template", Json{});

            //  ── primitives ─────────────────────────────
            if (typ == "int") {
                addPrimitive<int>(key,"int",defRaw,req,desc);
            }
            else if (typ == "bool") {
                addPrimitive<bool>(key,"bool",defRaw,req,desc);
            }
            else if (typ == "float") {
                addPrimitive<float>(key,"float",defRaw,req,desc);
            }
            else if (typ == "double") {
                addPrimitive<double>(key,"double",defRaw,req,desc);
            }
            else if (typ == "string") {
                addPrimitive<std::string>(key,"string",defRaw,req,desc);
            }
            //  ── registered Serializable types ───────
            else {
                auto& facs = _factories();
                if (!subTpl.is_null() || facs.count(typ)) {
                    // use custom Serializable factory
                    FieldInfo fi;
                    fi.name           = key;
                    fi.typeName       = typ;
                    fi.description    = desc;
                    fi.defaultRaw     = defRaw;
                    fi.required       = req;
                    fi.nestedTemplate = subTpl;

                    // create instance
                    {
                        std::lock_guard<std::mutex> lock(_reg_mtx());
                        if (!facs.count(typ))
                            throw std::runtime_error("Unregistered type: "+typ);
                        fi.storage    = _factories()[typ]();
                        fi.destructor = _destructors()[typ];
                        fi.toJson     = _toJsonF()[typ];
                        fi.fromJson   = _fromJsonF()[typ];
                    }
                    // apply default, if any
                    if (!defRaw.empty() && defRaw!="null") {
                        Json d = Json::parse(defRaw);
                        if (!fi.fromJson(fi.storage,d))
                            throw std::runtime_error("Bad default for "+key);
                    }
                    _fields[key] = std::move(fi);
                }
                else {
                    throw std::runtime_error("Unknown type: " + typ);
                }
            }
        }
        return true;
    }

    //───────────────────────────────────────────────────────────────────────────
    //  3) SerializableBase overrides
    //───────────────────────────────────────────────────────────────────────────
    Json GetTemplate() const override {
        Json out = Json::object();
        for (auto const& [k,fi] : _fields) {
            auto& o = out[k];
            o["type"]        = fi.typeName;
            o["description"] = fi.description;
            o["required"]    = fi.required;
            o["default"]     = fi.defaultRaw;
            if (!fi.nestedTemplate.is_null())
                o["template"] = fi.nestedTemplate;
        }
        return out;
    }
    std::string GetTemplateString() const override {
        return GetTemplate().dump();
    }

    Json GetJson() const override {
        Json out = Json::object();
        for (auto const& [k,fi] : _fields)
            out[k] = fi.toJson(fi.storage);
        return out;
    }
    std::string GetJsonString() const override {
        return GetJson().dump();
    }

    bool SetFromJson(const Json& j) override {
        for (auto const& [k,fi] : _fields) {
            auto it = j.find(k);
            if (it!=j.end()) {
                if (!fi.fromJson(fi.storage,*it))
                    return false;
            }
            else if (fi.required) {
                return false;
            }
        }
        return true;
    }
    bool SetFromJsonString(const std::string& s) override {
        return SetFromJson(Json::parse(s));
    }

    std::string GetTypeNameString() const override {
        return "DynamicSerializable";
    }

private:
    std::map<std::string,FieldInfo> _fields;

    void clear() {
        for (auto& [k,fi] : _fields) {
            if (fi.destructor) fi.destructor(fi.storage);
        }
        _fields.clear();
    }

    //───────────────────────────────────────────────────────────────────────────
    //  primitives helper
    template<typename T>
    void addPrimitive(const std::string& name,
                      const std::string& typeName,
                      const std::string& defaultRaw,
                      bool               required,
                      const std::string& description)
    {
        FieldInfo fi;
        fi.name        = name;
        fi.typeName    = typeName;
        fi.description = description;
        fi.defaultRaw  = defaultRaw;
        fi.required    = required;

        // native storage
        T* storage = new T();
        fi.storage    = storage;
        fi.destructor = [](void*p){ delete static_cast<T*>(p); };
        fi.toJson     = [](void*p){ return Json(*static_cast<T*>(p)); };
        fi.fromJson   = [](void*p,const Json& j){
            try {
                *static_cast<T*>(p) = j.get<T>();
                return true;
            } catch(...) {
                return false;
            }
        };

        // default
        if (!defaultRaw.empty() && defaultRaw!="null") {
            Json d = Json::parse(defaultRaw);
            if (!fi.fromJson(fi.storage,d))
                throw std::runtime_error("Bad default for "+name);
        }
        _fields[name] = std::move(fi);
    }

    //───────────────────────────────────────────────────────────────────────────
    //  static registry for Serializable types
    static std::map<std::string, std::function<void*()>>&
    _factories() {
        static std::map<std::string,std::function<void*()>> m;
        return m;
    }
    static std::map<std::string,std::function<void(void*)>>&
    _destructors() {
        static std::map<std::string,std::function<void(void*)>> m;
        return m;
    }
    static std::map<std::string,std::function<Json(void*)>>&
    _toJsonF() {
        static std::map<std::string,std::function<Json(void*)>> m;
        return m;
    }
    static std::map<std::string,std::function<bool(void*,const Json&)>>&
    _fromJsonF() {
        static std::map<std::string,std::function<bool(void*,const Json&)>> m;
        return m;
    }
    static std::mutex& _reg_mtx() {
        static std::mutex m;
        return m;
    }
};

#endif // DYNAMIC_SERIALIZABLE_HPP
