#ifndef CsfMapH
#define CsfMapH

#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <iostream>
#include <functional>
#include <atomic>

#define FF_OMP 1

namespace pcr {
    template<typename T>
    inline bool isMV(T val) { return val < -1e30f; }
    template<typename T>
    inline void setMV(T &val) { val = -1e31f; }
}

typedef float REAL4;
typedef unsigned char BYTE;

struct BoundingBox {
    double minx, maxx, miny, maxy;
    BoundingBox(double x1=0, double x2=0, double y1=0, double y2=0) : minx(x1), maxx(x2), miny(y1), maxy(y2) {}
};

struct GeoProjection {
    std::string wkt;
    static GeoProjection FromString(std::string s) { return GeoProjection{s}; }
    static GeoProjection GetFromEPSG(int epsg) { return GeoProjection{""}; }
    static GeoProjection GetFromWKT(std::string wkt) { return GeoProjection{wkt}; }
    std::string GetWKT() const { return wkt; }
    bool IsValid() const { return true; }
    bool IsProjected() const { return true; }
    int GetEPSG() const { return 3857; }
};

struct LSMVector2 {
    double x, y;
    LSMVector2(double x=0, double y=0) : x(x), y(y) {}
};

struct GeoCoordTransformer {
    static std::shared_ptr<GeoCoordTransformer> GetCoordTransformerSP(GeoProjection pfrom, GeoProjection pto) {
        return nullptr;
    }
    LSMVector2 Transform(LSMVector2 v) { return v; }
};

namespace LISEMMath {
    inline float CubicHermite(float A, float B, float C, float D, float t) {
        float a = -A / 2.0f + (3.0f * B) / 2.0f - (3.0f * C) / 2.0f + D / 2.0f;
        float b = A - (5.0f * B) / 2.0f + 2.0f * C - D / 2.0f;
        float c = -A / 2.0f + C / 2.0f;
        float d = B;
        return a * t * t * t + b * t * t + c * t + d;
    }
}

template<typename T>
class MaskedRaster {
public:
    std::vector<T> _data;
    size_t _nr_rows, _nr_cols;
    double _north, _west, _cell_sizeX, _cell_sizeY;

    MaskedRaster(size_t r=0, size_t c=0, double n=0, double w=0, double dx=1, double dy=1)
        : _nr_rows(r), _nr_cols(c), _north(n), _west(w), _cell_sizeX(dx), _cell_sizeY(dy) {
        _data.resize(r * c, 0);
    }
    size_t nr_rows() const { return _nr_rows; }
    size_t nr_cols() const { return _nr_cols; }
    size_t nr_cells() const { return _nr_rows * _nr_cols; }
    double north() const { return _north; }
    double west() const { return _west; }
    double cell_sizeX() const { return _cell_sizeX; }
    double cell_sizeY() const { return _cell_sizeY; }
    double cell_size() const { return _cell_sizeX; }

    T* operator[](size_t r) { return &_data[r * _nr_cols]; }
    const T* operator[](size_t r) const { return &_data[r * _nr_cols]; }

    T& cell(size_t index) { return _data[index]; }
    const T& cell(size_t index) const { return _data[index]; }
    T& cell(size_t r, size_t c) { return _data[r * _nr_cols + c]; }
    const T& cell(size_t r, size_t c) const { return _data[r * _nr_cols + c]; }

    void set_all_mv() {
        for(auto& v : _data) pcr::setMV(v);
    }
    void replace_with_mv(T const& value) {
        for(auto& v : _data) if(v == value) pcr::setMV(v);
    }
    void offsetscale(T const& scale, T const& offset) {
        for(auto& v : _data) if(!pcr::isMV(v)) v = v * scale + offset;
    }
};

template<typename T>
inline std::string to_string(T val) { return std::to_string(val); }

#define FOR_ROW_COL_MV(map) for(int r = 0; r < map->nrRows(); r++)\
    for (int c = 0; c < map->nrCols(); c++)\
    if(!pcr::isMV(map->data[r][c]))

#define FOR_ROW_COL(map) for(int r = 0; r < map->nrRows(); r++)\
    for (int c = 0; c < map->nrCols(); c++)\

#define MAPMATH_PI 3.14159f

#define Drc data[r][c]

#define LISEM_ASMAP_BOOLFROMFLOAT(x) (x>0.5f?true:false)

#define SWAPMAP(x, y) do { cTMap * SWAP = x; x = y; y = SWAP; } while (0)

static inline bool FuzzyCompare(double p1, double p2)
{
    return (std::abs(p1 - p2) * 1000000000000. <= std::min(std::abs(p1), std::abs(p2)));
}

static inline bool FuzzyCompare(float p1, float p2)
{
    return (std::abs(p1 - p2) * 100000.f <= std::min(std::abs(p1), std::abs(p2)));
}

//forward declaration of class
template<typename T>
class cTMap;

//forward declaration of map arithmatic base function
template<typename T, typename G, typename Callable>
cTMap<T> MapArithmatic(const cTMap<T>&,const cTMap<T>&,const Callable &F);
//forward declaration of map arithmatic base function
template<typename T, typename G, typename Callable>
cTMap<T> MapRealArithmatic(const cTMap<T>&,const G&,const Callable &F);

//forward declaration of map arithmatic base function
template<typename T, typename Callable>
void MapArithmaticInPlace(cTMap<T>&,const cTMap<T>&,const Callable &F);
//forward declaration of map arithmatic base function
template<typename T, typename Callable>
void MapRealArithmaticInPlace(cTMap<T>&,const T&,const Callable &F);

template<typename T,  typename Callable>
cTMap<T> CallMapSingleFunc(cTMap<T> * m,const Callable &f);

typedef struct RasterBandStats
{
    double min = 0.0;
    double max = 0.0;
    double stdev = 0.0;
    double mean = 0.0;
    int band = 0;
    int n_mv = 0;
    int rows = 1;
    int cols = 1;
    double dx = 1.0;
    double dy = 1.0;
} RasterBandStats;

typedef struct cTMapProps
{
    bool can_be_read = false;
    std::string file_path = "";
    int sizex = 0.0;
    int sizey = 0.0;
    int bands = 0;
    int band = 0;
    bool is_ldd = false;
    double cellsizex = 0.0;
    double cellsizey = 0.0;

    std::string projection;

    double uly = 0.0;
    double ulx = 0.0;

    bool equals(cTMapProps p);

} cTMapProps;


/*!
    @brief      A cTMap contains all relevant information about a raster.
*/

template<typename T = float>
class cTMap
{

public:

std::string m_projection = "";


    //! The actual raster.
    MaskedRaster<T> data;
    BoundingBox m_BoundingBox;

    std::string m_FileName =  "";
    bool m_IsLDD = false;
                   cTMap               ()=default;

                   cTMap 			   (float * datain,
                                        int width, int height,
                                        float north, float west,
                                        float dx, float dy,
                                        std::string projection ="");

                   cTMap               (MaskedRaster<float>&& data,
                                        bool is_ldd = false);

                   cTMap               (MaskedRaster<float>&& data,
										std::string projection,
                                        std::string mapName,
                                        bool is_ldd = false);

                   cTMap               (cTMap<T> const& other)=default;

                   cTMap               (cTMap<T> && other)=default;

                   ~cTMap              (){}//=default;

    inline cTMap&  operator=           (cTMap const& other)
    {
        //assignment operator,
        //what if this object is reference counted,
        //it should not copy the reference count as well
        data = other.data;
        m_IsLDD = other.m_IsLDD;
        m_FileName = other.m_FileName;
        return *this;
    }

    cTMap&         operator=           (cTMap&& other)=default;

    int            nrRows              () const;

    int            nrCols              () const;

    double         north               () const;
    double         west                () const;
    double         cellSize            () const;
    double         cellSizeX            () const;
    double         cellSizeY            () const;

    double CellSizeXApproximateMeters();
    void           setnorth            (double);
    void           setwest             (double);
    void           setcellSizeX        (double);
    void           setcellSizeY        (double);

    T              val_min             ();
    T              val_max             ();

    T              DerX                (int r, int c);
    T              DerY                (int r, int c);
    void           setAllMV            ();
    void           setAllValue         (T value);

    void           MakeMap             (cTMap *dup,
                                        REAL4 value);

    LSMVector2     GetPixelLngLat       (int r, int c);

    T              GetPixelCovered      (int r, int c, float cover = 0.0);
    T              Sample              (double x, double y);
    T              SampleLinear              (double x, double y);
    T              SampleLinearPS              (double x, double y);
    T              SampleCubic              (double x, double y);
    cTMap*         GetCopy             () const;
    cTMap*         GetCopy0             () const;
    std::shared_ptr<cTMap>    GetCopySP0             () const;
    std::shared_ptr<cTMap>    GetCopySP             () const;
    cTMap*         GetCopy0Cascade             (int pixel_scale_mult) const;
    cTMap*         GetCopy1             () const;
    cTMap*         GetCopy1Masked             () const;
    cTMap*         GetCopyMask          () const;

    bool                CheckIsEdgeOfMaskedDomain   (size_t r, size_t c);
    std::string         projectionstring                  ()                      const;
    void                setprojectionstring               (std::string proj);
    BoundingBox         GetExtent                   ();
    RasterBandStats     GetRasterBandStats          (bool stdev = false);

    inline void SetProjection(std::string p) { m_projection = p; }

    //iterators and loops
public:

    inline virtual GeoProjection GetProjection()
    {
        return GeoProjection::FromString(m_projection);
    }

    template<typename Callable>
    inline void ForEach(Callable &&f)
    {
        cTMap<T> * map1=  this;
        #if(FF_OMP)
            #pragma omp parallel for collapse(2)
            for(int r = 0; r < map1->data.nr_rows();r++)
            {
                for(int c = 0; c < map1->data.nr_cols();c++)
                {
                    if(!pcr::isMV(map1->data[r][c]))
                    {
                        f(r,c);
                    }
                }
            }

        #else
            GetTPool().parallel_for(0,map1->data.nr_rows(),[&](int r)
            {
                for(int c = 0; c < map1->data.nr_cols();c++)
                {
                    if(!pcr::isMV(map1.data[r][c]))
                    {
                        f(r,c);
                    }
                }

            });

        #endif

    }
    //Operators for arithmatic
public:

    inline cTMap<T> &Assign(const cTMap &in)
    {
        *this= in;
        return *this;
    }

    //template<typename G>
    inline cTMap<T> operator+(const cTMap<T> & in) const
    {
        return MapArithmatic(*this,in,[](const T & l, const T & r){
            return (T)(l+r);
        });
    }
    inline cTMap<T>* OpAdd(const cTMap<T> * in) const
    {
        cTMap<T> ret =MapArithmatic(*this,*in,[](const T & l, const T & r)->T{
            return (T)(l+r);
        });
        return new cTMap<T>(std::move(ret));
    }

    inline cTMap<T> operator-(const cTMap<T> & in) const
    {
        return MapArithmatic(*this,in,[](const T & l, const T & r){
            return (T)(l-r);
        });
    }
    inline cTMap<T>* OpSub(const cTMap<T> * in) const
    {
        return new cTMap<T>(std::move(MapArithmatic(*this,*in,[](const T & l, const T & r){
            return (T)(l-r);
        })));
    }

    inline cTMap<T> operator*(const cTMap<T> & in) const
    {
        return MapArithmatic(*this,in,[](const T & l, const T & r){
            return (T)(l+r);
        });
    }
    inline cTMap<T>* OpMul(const cTMap<T> * in) const
    {
        return new cTMap<T>(std::move(MapArithmatic(*this,*in,[](const T & l, const T & r){
            return (T)(l+r);
        })));
    }
    inline cTMap<T> operator/(const cTMap<T> & in) const
    {
        return MapArithmatic(*this,in,[](const T & l, const T & r){
            return (T)(l/r);
        });
    }
    inline cTMap<T>* OpDiv(const cTMap<T> * in) const
    {
        return new cTMap<T>(std::move(MapArithmatic(*this,*in,[](const T & l, const T & r){
            return (T)(l/r);
        })));
    }

    inline cTMap<T> operator^(const cTMap<T> & in) const
    {
        return MapArithmatic(*this,in,[](const T & l, const T & r){
            return (T)(std::pow(l,r));
        });
    }
    inline cTMap<T>* OpPow(const cTMap<T> * in) const
    {
        return new cTMap<T>(std::move(MapArithmatic(*this,*in,[](const T & l, const T & r){
            return (T)(std::pow(l,r));
        })));
    }

    inline cTMap<T> operator%(const cTMap<T> & in) const
    {
        return MapArithmatic(*this,in,[](const T & l, const T & r){
            return (T)(std::fmod(l,r));
        });
    }
    inline cTMap<T>* OpMod(const cTMap<T> * in) const
    {
        return new cTMap<T>(std::move(MapArithmatic(*this,*in,[](const T & l, const T & r){
            return (T)(std::fmod(l,r));
        })));
    }

    inline cTMap<T> operator+(const T & in) const
    {
        return MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(l+r);
        });
    }
    inline cTMap<T>* OpAdd(const T & in) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(l+r);
        })));
    }

    inline cTMap<T> operator-(const T & in) const
    {
        return MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(l-r);
        });
    }
    inline cTMap<T>* OpSub(const T & in) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(l-r);
        })));
    }

    inline cTMap<T> operator*(const T & in) const
    {
        return MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(l*r);
        });
    }
    inline cTMap<T>* OpMul(const T & in) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(l*r);
        })));
    }
    inline cTMap<T> operator/(const T & in) const
    {
        return MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(l*r);
        });
    }
    inline cTMap<T>* OpDiv(const T & in) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(l/r);
        })));
    }
    inline cTMap<T> operator^(const T & in) const
    {
        return MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(std::pow(l,r));
        });
    }
    inline cTMap<T>* OpPow(const T & in) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(std::powf(l,r));
        })));
    }
    inline cTMap<T> operator%(const T & in) const
    {
        return MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(std::fmod(l,r));
        });
    }
    inline cTMap<T>* OpMod(const T & in) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic<T,T>(*this,in,[](const T & l, const T & r){
            return (T)(std::fmod(l,r));
        })));
    }


    //template<typename G>
    inline cTMap<T> &operator+=(const cTMap<T> & in)
    {
        MapArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l+r);
        });
        return this;
    }
    inline cTMap<T>* OpAddAssign(const cTMap<T> * in)
    {
        MapArithmaticInPlace(*this,*in,[](const T & l, const T & r){
            return (T)(l+r);
        });
        return this;
    }

    inline cTMap<T> &operator-=(const cTMap<T> & in)
    {
        MapArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l-r);
        });
        return *this;
    }
    inline cTMap<T>* OpSubAssign(const cTMap<T> * in)
    {
        MapArithmaticInPlace(*this,*in,[](const T & l, const T & r){
            return (T)(l-r);
        });
        return this;
    }

    inline cTMap<T> &operator*=(const cTMap<T> & in)
    {
        MapArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l*r);
        });
        return *this;
    }
    inline cTMap<T>* OpMulAssign(const cTMap<T> * in)
    {
        MapArithmaticInPlace(*this,*in,[](const T & l, const T & r){
            return (T)(l*r);
        });
        return this;
    }
    inline cTMap<T> &operator/=(const cTMap<T> & in)
    {
        MapArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l/r);
        });
        return *this;
    }
    inline cTMap<T>* OpDivAssign(const cTMap<T> * in)
    {
        MapArithmaticInPlace(*this,*in,[](const T & l, const T & r){
            return (T)(l/r);
        });
        return this;
    }

    inline cTMap<T> &operator^=(const cTMap<T> & in)
    {
        MapArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(std::pow(l,r));
        });
        return *this;
    }
    inline cTMap<T>* OpPowAssign(const cTMap<T> * in)
    {
        MapArithmaticInPlace(*this,*in,[](const T & l, const T & r){
            return (T)(std::pow(l,r));
        });
        return this;
    }

    inline cTMap<T> &operator%=(const cTMap<T> & in)
    {
        MapArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(std::fmod(l,r));
        });
        return *this;
    }
    inline cTMap<T>* OpModAssign(const cTMap<T> * in)
    {
        MapArithmaticInPlace(*this,*in,[](const T & l, const T & r){
            return (T)(std::fmod(l,r));
        });
        return this;
    }

    inline cTMap<T> &operator+=(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l+r);
        });
        return *this;
    }
    inline cTMap<T>* OpAddAssign(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l+r);
        });
        return this;
    }

    inline cTMap<T> &operator-=(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l-r);
        });
        return *this;
    }
    inline cTMap<T>* OpSubAssign(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l-r);
        });
        return this;
    }

    inline cTMap<T> &operator*=(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l*r);
        });
        return *this;
    }
    inline cTMap<T>* OpMulAssign(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l*r);
        });
        return this;
    }
    inline cTMap<T> &operator/=(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l/r);
        });
        return *this;
    }
    inline cTMap<T>* OpDivAssign(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(l/r);
        });
        return this;
    }
    inline cTMap<T> &operator^=(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(std::powf(l,r));
        });
        return *this;
    }
    inline cTMap<T>* OpPowAssign(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(std::powf(l,r));
        });
        return this;
    }
    inline cTMap<T> &operator%=(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(std::fmod(l,r));
        });
        return *this;
    }
    inline cTMap<T>* OpModAssign(const T & in)
    {
        MapRealArithmaticInPlace(*this,in,[](const T & l, const T & r){
            return (T)(std::fmod(l,r));
        });
        return this;
    }

    inline cTMap<T> * EqualTo(const cTMap<T> *other) const
    {
        std::function<T(T,T)> eqfunc = [](T a, T b){return (FuzzyCompare(a,b)? 1.0f: 0.0f);};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,eqfunc)));
    }
    inline cTMap<T> * LargerThen(const cTMap<T> *other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a > b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,compfunc)));
    }
    inline cTMap<T> * SmallerThen(const cTMap<T> *other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a < b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,compfunc)));
    }
    inline cTMap<T> * NotEqualTo(const cTMap<T> *other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a != b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,compfunc)));
    }
    inline cTMap<T> * LargerEqualThen(const cTMap<T> *other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a >= b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,compfunc)));
    }
    inline cTMap<T> * SmallerEqualThen(const cTMap<T> *other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a <= b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,compfunc)));
    }
    inline cTMap<T> * EqualTo(const T &other) const
    {

        std::function<T(T,T)> compfunc = [](T a, T b){return (FuzzyCompare(a,b)? 1.0f: 0.0f);};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * LargerThen(const T &other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a > b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));

    }
    inline cTMap<T> * SmallerThen(const T &other) const
    {

        std::function<T(T,T)> compfunc = [](T a, T b){return (a < b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * NotEqualTo(const T &other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a != b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * LargerEqualThen(const T &other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a >= b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * SmallerEqualThen  (const T &other) const
    {
        std::function<T(T,T)> compfunc = [](T a, T b){return (a <= b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }

    inline cTMap<T> * And               (const T &other) const
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) && LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }
    inline cTMap<T> * Or                (const T &other) const
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) || LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }
    inline cTMap<T> * Xor               (const T &other) const
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) ^ LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }

    inline cTMap<T> * And               (const cTMap * other) const
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) && LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,func)));
    }
    inline cTMap<T> * Or                (const cTMap * other) const
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) || LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,func)));

    }
    inline cTMap<T> * Xor               (const cTMap * other) const
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) ^ LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapArithmatic(*this,*other,func)));
    }
    inline cTMap<T> * AndAssign              (const T &other)
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) && LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        MapRealArithmaticInPlace(*this,other,func);
        return this;
    }
    inline cTMap<T> * OrAssign                (const T &other)
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) || LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        MapRealArithmaticInPlace(*this,other,func);
        return this;
    }
    inline cTMap<T> * XorAssign               (const T &other)
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) ^ LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        MapRealArithmaticInPlace(*this,other,func);
        return this;
    }
    inline cTMap<T> * AndAssign               (const cTMap * other)
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) && LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        MapArithmaticInPlace(*this,*other,func);
        return this;
    }
    inline cTMap<T> * OrAssign                (const cTMap * other)
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) || LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        MapArithmaticInPlace(*this,*other,func);
        return this;

    }
    inline cTMap<T> * XorAssign               (const cTMap * other)
    {
        std::function<T(T,T)> func = [](T a, T b){return (LISEM_ASMAP_BOOLFROMFLOAT(a) ^ LISEM_ASMAP_BOOLFROMFLOAT(b))? 1.0f:0.0f;};
        MapArithmaticInPlace(*this,*other,func);
        return this;
    }
    inline cTMap<T> * Negate () const
    {
        std::function<T(T,T)> func = [](T a, T b){return (!(LISEM_ASMAP_BOOLFROMFLOAT(a)) )? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,1,func)));
    }
    inline cTMap<T> * OpMod_r(const T& other) const
    {
        std::function<T(T,T)> func = [](T a, T b){return std::fmod(b,a);};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }

    inline cTMap<T> * OpAdd_r(const T& other) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,std::plus<T>())));
    }
    inline cTMap<T> * OpMul_r(const T& other) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,std::multiplies<T>())));
    }
    inline cTMap<T> * OpSub_r(const T& other) const
    {
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,[](T a, T b){return b-a;})));
    }
    inline cTMap<T> * OpDiv_r(const T& other) const
    {
        std::function<T(T,T)> func = [](const T& a, const T& b){return b/a;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }
    inline cTMap<T> * OpPow_r(const T& other) const
    {
        std::function<T(T,T)> func = [](const T& a, const T& b){return std::pow(b,a);};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }
    inline cTMap<T> * LargerThen_r      (const T &other) const
    {
        std::function<T(T,T)> compfunc = [](const T& a, const T& b){return (b > a)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * SmallerThen_r     (const T &other) const
    {
        std::function<T(T,T)> compfunc = [](const T& a, const T& b){return (b < a)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * NotEqualTo_r      (const T& other) const
    {
        std::function<T(T,T)> compfunc = [](const T& a, const T& b){return (a != b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * LargerEqualThen_r (const T& other) const
    {
        std::function<T(T,T)> compfunc = [](const T& a, const T& b){return (a >= b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * SmallerEqualThen_r(const T& other) const
    {
        std::function<T(T,T)> compfunc = [](const T& a, const T& b){return (a <= b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * EqualTo_r         (const T& other) const
    {
        std::function<T(T,T)> compfunc = [](const T& a, const T& b){return (a == b)? 1.0f: 0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,compfunc)));
    }
    inline cTMap<T> * And_r (const T& other) const
    {
        std::function<T(T,T)> func = [](const T& a, const T& b){return (LISEM_ASMAP_BOOLFROMFLOAT(b) && LISEM_ASMAP_BOOLFROMFLOAT(a))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }
    inline cTMap<T> * Or_r  (const T& other) const
    {
        std::function<T(T,T)> func = [](const T& a, const T& b){return (LISEM_ASMAP_BOOLFROMFLOAT(b) || LISEM_ASMAP_BOOLFROMFLOAT(a))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }
    inline cTMap<T> * Xor_r (const T& other) const
    {
        std::function<T(T,T)> func = [](const T& a, const T& b){return (LISEM_ASMAP_BOOLFROMFLOAT(b) ^ LISEM_ASMAP_BOOLFROMFLOAT(a))? 1.0f:0.0f;};
        return new cTMap<T>(std::move(MapRealArithmatic(*this,other,func)));
    }

    inline virtual void OnAssign()
    {

    }
public:

    float&         Value (int r, int c);
    int            Rows              ();
    int            Cols              ();
    float          North             ();
    float          West              ();
    float          CellSize          ();
    cTMapProps     GetRasterProps      ();



};

//typedefs for commonly used template types

typedef cTMap<float> cTMapf;
typedef cTMap<double> cTMapd;
typedef cTMap<int> cTMapi;
typedef cTMap<long int> cTMapl;
typedef cTMap<unsigned long int> cTMapul;
typedef cTMap<unsigned int> cTMapui;


//serialization
template<typename T>
inline std::string to_string(std::shared_ptr<cTMap<T>> data, bool do_compress = true) { return ""; }
inline std::shared_ptr<cTMapf> from_string(std::string x) { return nullptr; }

//some base global functions to deal with cTMap manipulation, non-algebraic or algorithmic related
template<typename T>
static void DeNormalizeMapInPlace(cTMap<T> *m, double min, double max);
template<typename T>
static void NormalizeMapInPlace(cTMap<T> *m, double &min, double &max);

template<typename T>
static cTMap<T> * NormalizeMap(cTMap<T> *m, double &min, double &max);

template<typename T>
inline bool OUTORMV(cTMap<T> * map, int r, int c);



//==========================================================
//                  Implementation details for
//                  templated map class
//==========================================================

template<typename T>
inline bool cTMapProps::equals(cTMapProps p)
{
    if(p.can_be_read != can_be_read)
    {
        return false;
    }
    if(sizex != p.sizex)
    {
        return false;
    }
    if(sizey != p.sizey)
    {
        return false;
    }
    if(bands != p.bands)
    {
        return false;
    }
    if(is_ldd != p.is_ldd)
    {
        return false;
    }
    if(cellsizex != p.cellsizex)
    {
        return false;
    }
    if(cellsizey != p.cellsizey)
    {
        return false;
    }
    if(projection != p.projection)
    {
        return false;
    }
    if(ulx != p.ulx)
    {
        return false;
    }
    if(uly != p.uly)
    {
        return false;
    }
    return true;
}


template<typename T>
inline cTMapProps     cTMap<T>::GetRasterProps      ()
{
    cTMapProps p;

    p.ulx = west();
    p.uly = north();
    p.band = 0;
    p.bands = 1;
    p.sizex = nrCols();
    p.sizey = nrRows();
    p.is_ldd = m_IsLDD;
    p.projection = projectionstring();
    p.cellsizex = cellSizeX();
    p.cellsizey = cellSizeY();
    p.can_be_read = m_FileName.size() == 0?false:true;
    p.file_path = m_FileName;

    return p;
}

template<typename T>
inline RasterBandStats cTMap<T>::GetRasterBandStats(bool do_stdev)
{
    RasterBandStats r;

    float hmax = -1e31f;
    float hmin = 1e31f;
    float havg = 0.0f;
    float n = 0.0;
    float stdev = 0.0;

    int n_mv = 0;
    for(int row = 0; row < this->nrRows(); row++)
    {
        for (int col = 0; col < this->nrCols(); col++)
        {


        if(!pcr::isMV(this->data[row][col]))
        {
            if(std::isfinite(this->data[row][col]) && !std::isnan((this->data[row][col])))
            {
                n++;
                havg += (this->data[row][col]);
                hmax = std::max(hmax,this->data[row][col]);
                hmin = std::min(hmin,this->data[row][col]);
            }else {
                n_mv ++;
            }
        }else {
            n_mv++;
        }
        }
    }

    n = std::max(1.0f,n);

    if(do_stdev)
    {
        FOR_ROW_COL_MV(this)
        {
            if(std::isfinite(this->data[r][c]) && !std::isnan((this->data[r][c])))
            {
                stdev += (this->data[r][c] -havg)*(this->data[r][c] -havg);
            }
        }
        stdev = stdev/n;
    }else {
        stdev = 0.35 * (hmax-hmin);
    }

    r.mean = havg/n;
    r.min = hmin;
    r.max = hmax;
    r.stdev = std::sqrt(stdev);
    r.n_mv = n_mv;

    return r;

}
template<typename T>
inline cTMap<T>::cTMap(float * datain, int width, int height, float innorth, float inwest, float dx,float dy, std::string projection)
{
    data = MaskedRaster<float>(height,width, innorth,inwest,dx,dy);
    if(datain != nullptr)
    {

        for(int r = 0; r < height; r++)
        {
            for(int c = 0; c < width; c++)
            {
                int index = r * width + c;
                if(datain[index] < -1e30)
                {
                    pcr::setMV(data[r][c]);
                }else{

                data[r][c] = datain[index];
                }
            }
        }
    }else{
        for(int r = 0; r < height; r++)
        {
            for(int c = 0; c < width; c++)
            {

                data[r][c] = 0.0;

            }
        }

    }

    this->m_projection = projection;
    this->SetProjection(m_projection);
    m_BoundingBox = BoundingBox(west(),west() + cellSizeX()*((double)(nrCols())),north(),north() + cellSizeY()*((double)(nrRows())));
}

template<typename T>
inline cTMap<T>::cTMap(
    MaskedRaster<float>&& data,
        bool is_ldd)

    : data(std::forward<MaskedRaster<float>>(data))

{

    m_BoundingBox = BoundingBox(west(),west() + cellSizeX()*((double)(nrCols())),north(),north() + cellSizeY()*((double)(nrRows())));
}
template<typename T>
inline cTMap<T>::cTMap(
    MaskedRaster<float>&& data,
    std::string project,
    std::string name,
        bool is_ldd)

    : data(std::forward<MaskedRaster<float>>(data))

{

    m_projection = project;
    SetProjection(m_projection);
    m_BoundingBox = BoundingBox(west(),west() + cellSizeX()*((double)(nrCols())),north(),north() + cellSizeY()*((double)(nrRows())));

}
template<typename T>
inline int cTMap<T>::nrRows() const
{
    return static_cast<int>(data.nr_rows());
}

template<typename T>
inline int cTMap<T>::nrCols() const
{
    return static_cast<int>(data.nr_cols());
}
template<typename T>
inline double cTMap<T>::north() const
{
    return data.north();
}
template<typename T>
inline double cTMap<T>::west() const
{
    return data.west();
}

template<typename T>
inline double cTMap<T>::cellSize() const
{
    return data.cell_size();
}
template<typename T>
inline double cTMap<T>::cellSizeX() const
{
    return data.cell_sizeX();
}

template<typename T>
inline double cTMap<T>::CellSizeXApproximateMeters()
{
    double dx = std::fabs(data.cell_sizeX());
    GeoProjection p = GeoProjection::GetFromWKT(projectionstring());

    if(p.IsValid() && !p.IsProjected())
    {
        float ymin = North();
        float ymax = North() + (((float)Rows()) * cellSizeY());
        dx = dx * 111321.0 * std::cos(0.5* (ymin + ymax) * 3.14159/180.0);
    }

    return dx;
}

template<typename T>
inline double cTMap<T>::cellSizeY() const
{
    return data.cell_sizeY();
}
template<typename T>
inline void cTMap<T>::setnorth(double x)
{
    data._north = x;

}
template<typename T>
inline void cTMap<T>::setwest(double x)
{

    data._west = x;
}
template<typename T>
inline void cTMap<T>::setcellSizeX(double x)
{

    data._cell_sizeX= x;
    m_BoundingBox = BoundingBox(west(),west() + cellSizeX()*((double)(nrCols())),north(),north() + cellSizeY()*((double)(nrRows())));
}
template<typename T>
inline void cTMap<T>::setcellSizeY(double x)
{
    data._cell_sizeY= x;
    m_BoundingBox = BoundingBox(west(),west() + cellSizeX()*((double)(nrCols())),north(),north() + cellSizeY()*((double)(nrRows())));
}
template<typename T>
inline void cTMap<T>::setAllMV()
{
    data.set_all_mv();
}
template<typename T>
inline void cTMap<T>::setAllValue(T val)
{
    this->ForEach([this,val](size_t r, size_t c)
    {
        this->data[r][c] = val;
    });
}

// make a new map according to dup as a mask and filled with value
template<typename T>
inline void cTMap<T>::MakeMap(
    cTMap *dup,
    REAL4 value)
{
  if (dup == nullptr)
    return;

  data = MaskedRaster<REAL4>(dup->nrRows(), dup->nrCols(), dup->north(),
      dup->west(), dup->cellSizeX(),dup->cellSizeY());

  data.set_all_mv();

  for(int r=0; r < nrRows(); r++)
    for(int c=0; c < nrCols(); c++)
      if (!pcr::isMV(dup->data[r][c]))
        {
          data[r][c] = value;
        }
  m_BoundingBox = BoundingBox(west(),west() + cellSizeX()*((double)(nrCols())),north(),north() + cellSizeY()*((double)(nrRows())));

}
template<typename T>
inline cTMap<T> * cTMap<T>::GetCopy() const
{
    const cTMap<T> * dup = this;

    cTMap<T> * m = new cTMap<T>(nullptr, dup->nrCols(),dup->nrRows(), dup->north(),
        dup->west(), dup->cellSizeX(),dup->cellSizeY(),this->projectionstring());

    m->setAllMV();

    for(int r=0; r < m->nrRows(); r++)
    {
      for(int c=0; c < m->nrCols(); c++)
      {
        if (!pcr::isMV(data[r][c]))
          {
            m->data[r][c] = data[r][c];
          }
      }
    }

    return m;
}
template<typename T>
inline cTMap<T> * cTMap<T>::GetCopyMask() const
{
    const cTMap<T> * dup = this;
    MaskedRaster<REAL4> datac = MaskedRaster<REAL4>(dup->nrRows(), dup->nrCols(), dup->north(),
        dup->west(), dup->cellSizeX(),dup->cellSizeY());

    cTMap<T> * m = new cTMap<T>(std::move(datac),false);
    m->setAllMV();

    for(int r=0; r < m->nrRows(); r++)
    {
      for(int c=0; c < m->nrCols(); c++)
      {
        if (!pcr::isMV(data[r][c]))
          {
            m->data[r][c] = 1.0;
          }else {
            m->data[r][c] = 0.0;
        }
      }
    }

    return m;
}
template<typename T>
inline std::shared_ptr<cTMap<T>> cTMap<T>::GetCopySP() const
{
    const cTMap<T> * dup = this;

    std::shared_ptr<cTMap<T>> sp = std::shared_ptr<cTMap<T>>(new cTMap<T>(nullptr, dup->nrCols(),dup->nrRows(), dup->north(),
        dup->west(), dup->cellSizeX(),dup->cellSizeY(),dup->projectionstring()));
    for(int r=0; r < sp->nrRows(); r++)
        {
          for(int c=0; c < sp->nrCols(); c++)
          {
            if (!pcr::isMV(data[r][c]))
              {
                sp->data[r][c] = data[r][c];
              }
          }
        }
    return sp;
}

template<typename T>
inline std::shared_ptr<cTMap<T>> cTMap<T>::GetCopySP0() const
{
    const cTMap<T> * dup = this;

    std::shared_ptr<cTMap<T>> sp = std::shared_ptr<cTMap<T>>(new cTMap<T>(nullptr, dup->nrCols(),dup->nrRows(), dup->north(),
        dup->west(), dup->cellSizeX(),dup->cellSizeY(),dup->projectionstring()));

    return sp;
}
template<typename T>
inline cTMap<T> * cTMap<T>::GetCopy0() const
{
    const cTMap<T> * dup = this;

    cTMap<T> * m = new cTMap<T>(nullptr, dup->nrCols(),dup->nrRows(), dup->north(),
        dup->west(), dup->cellSizeX(),dup->cellSizeY(),dup->projectionstring());

    return m;
}
template<typename T>
inline cTMap<T> * cTMap<T>::GetCopy0Cascade(int pixel_scale_mult) const
{
    const cTMap<T> * dup = this;
    int n_cols = dup->nrCols()/pixel_scale_mult;
    int n_rows = dup->nrRows()/pixel_scale_mult;
    if(n_cols * pixel_scale_mult < dup->nrCols())
    {
        n_cols +=1;
    }
    if(n_rows * pixel_scale_mult < dup->nrRows())
    {
        n_rows +=1;
    }
    cTMap<T> * m = new cTMap<T>(nullptr, n_cols,n_rows, dup->north(),
        dup->west(), dup->cellSizeX()*((float)pixel_scale_mult),dup->cellSizeY()*((float)pixel_scale_mult),this->projectionstring());

    return m;
}
template<typename T>
inline cTMap<T> * cTMap<T>::GetCopy1Masked() const
{
    const cTMap<T> * dup = this;

    cTMap<T> * m = new cTMap<T>(nullptr, dup->nrCols(),dup->nrRows(), dup->north(),
        dup->west(), dup->cellSizeX(),dup->cellSizeY(),this->projectionstring());


    for(int r=0; r < m->nrRows(); r++)
    {
      for(int c=0; c < m->nrCols(); c++)
      {
        if (!pcr::isMV(data[r][c]))
          {
            m->data[r][c] = 1.0;
          }else {
            pcr::setMV(m->data[r][c]);
        }
      }
    }

    return m;
}

template<typename T>
inline cTMap<T> * cTMap<T>::GetCopy1() const
{
    const cTMap<T> * dup = this;

    cTMap<T> * m = new cTMap<T>(nullptr, dup->nrCols(),dup->nrRows(), dup->north(),
        dup->west(), dup->cellSizeX(),dup->cellSizeY(),this->projectionstring());


    for(int r=0; r < m->nrRows(); r++)
    {
      for(int c=0; c < m->nrCols(); c++)
      {
        if (!pcr::isMV(data[r][c]))
          {
            m->data[r][c] = 1.0;
          }else {
            m->data[r][c] = 0.0;
        }
      }
    }

    return m;
}

template<typename T>
inline float& cTMap<T>::Value(int r, int c)
{
    return this->data[r][c];
}

template<typename T>
inline int cTMap<T>::Rows()
{
    return this->nrRows();
}
template<typename T>
inline int cTMap<T>::Cols()
{
    return this->nrCols();
}
template<typename T>
inline float cTMap<T>::North()
{
    return this->north();

}
template<typename T>
inline float cTMap<T>::West()
{
    return this->west();

}
template<typename T>
inline float cTMap<T>::CellSize()
{
    return this->cellSize();

}

template<typename T>
inline T cTMap<T>::SampleLinearPS(double x, double y)
{
    float ans = 0.0;

    float x_pix = x;
    float y_pix = y;

    //get closest pixel


    int xpix = x_pix;
    int ypix = y_pix;

    float x_frac = 1.0-(x_pix - (double)(xpix));
    float y_frac =  (cellSizeY() < 0)? 1.0-(y_pix - (double)(ypix)) : (y_pix - (double)(ypix));

    ypix = (cellSizeY() < 0? ypix + 1 : ypix);

    int xpix_x = xpix+(cellSizeX() > 0?1:-1);
    int ypix_y = ypix+(cellSizeY() > 0?1:-1);

    T val_x1y1,val_x1y2,val_x2y1,val_x2y2;

    int tmv =0;
    T tvaltemp = 0.0;
    int bmv =0;
    T bvaltemp = 0.0;

    if(!((xpix < 0)|| (xpix > (nrCols() - 1)) || (ypix < 0)|| (ypix > (nrRows() - 1))))
    {
       val_x1y1 = data[ypix][xpix];
       if(pcr::isMV(val_x1y1))
       {
           tmv ++;
       }else {
           tvaltemp = val_x1y1;
        }
    }else {
       pcr::setMV(val_x1y1);
       tmv ++;
    }

    if(!((xpix_x < 0)|| (xpix_x > (nrCols() - 1)) || (ypix < 0)|| (ypix > (nrRows() - 1))))
    {
       val_x2y1 = data[ypix][xpix_x];
       if(pcr::isMV(val_x2y1))
       {
           tmv ++;
       }else {
           tvaltemp = val_x2y1;
        }
    }else {
       pcr::setMV(val_x2y1);
       tmv ++;
    }

    if(!((xpix < 0)|| (xpix > (nrCols() - 1)) || (ypix_y < 0)|| (ypix_y > (nrRows() - 1))))
    {
       val_x1y2 = data[ypix_y][xpix];
       if(pcr::isMV(val_x1y2))
       {
           bmv ++;
       }else {
           bvaltemp = val_x1y2;
        }
    }else {
       pcr::setMV(val_x1y2);
       bmv ++;
    }

    if(!((xpix_x < 0)|| (xpix_x > (nrCols() - 1)) || (ypix_y < 0)|| (ypix_y > (nrRows() - 1))))
    {
       val_x2y2 = data[ypix_y][xpix_x];
       if(pcr::isMV(val_x2y2))
       {
           bmv ++;
       }else {
           bvaltemp = val_x2y2;
        }
    }else {
       pcr::setMV(val_x2y2);
       bmv++;
    }

    if((tmv == 2) && (bmv == 2))
    {
        float val;
        pcr::setMV(val);
        return val;

    }else if((tmv == 1) && (bmv == 2))
    {
        return tvaltemp;
    }else if((tmv == 2) && (bmv == 1))
    {
        return bvaltemp;
    }else if((tmv == 0) && (bmv == 2))
    {
        return x_frac * val_x2y1 + (1.0-x_frac)* val_x1y1;
    }else if((tmv == 2) && (bmv == 0))
    {
        return x_frac * val_x1y2 + (1.0-x_frac)* val_x2y2;
    }else if((tmv == 1) && (bmv == 0))
    {
        return y_frac * ( x_frac * val_x1y2 + (1.0-x_frac)* val_x2y2) + (1.0-y_frac)* tvaltemp;
    }else if((tmv == 1) && (bmv == 1))
    {
        return y_frac * bvaltemp + (1.0-y_frac)* tvaltemp;
    }else if((tmv == 0) && (bmv == 1))
    {
        return y_frac * bvaltemp + (1.0-y_frac)* ( x_frac * val_x1y1 + (1.0-x_frac)* val_x2y1) ;
    }else if((tmv == 0) && (bmv == 0))
    {
        return y_frac * ( x_frac * val_x1y2 + (1.0-x_frac)* val_x2y2) + (1.0-y_frac)* ( x_frac * val_x1y1 + (1.0-x_frac)* val_x2y1) ;
    }

    return data[ypix][xpix];
}
template<typename T>
inline T cTMap<T>::SampleLinear(double x, double y)
{
    float ans = 0.0;

    float x_pix = (x - west())/cellSizeX();
    float y_pix = (y - north())/cellSizeY();

    //get closest pixel


    int xpix = x_pix;
    int ypix = y_pix;

    float x_frac = 1.0-(x_pix - (double)(xpix));
    float y_frac =  (cellSizeY() < 0)? 1.0-(y_pix - (double)(ypix)) : (y_pix - (double)(ypix));

    ypix = (cellSizeY() < 0? ypix + 1 : ypix);

    int xpix_x = xpix+(cellSizeX() > 0?1:-1);
    int ypix_y = ypix+(cellSizeY() > 0?1:-1);

    T val_x1y1,val_x1y2,val_x2y1,val_x2y2;

    int tmv =0;
    T tvaltemp = 0.0;
    int bmv =0;
    T bvaltemp = 0.0;

    if(!((xpix < 0)|| (xpix > (nrCols() - 1)) || (ypix < 0)|| (ypix > (nrRows() - 1))))
    {
       val_x1y1 = data[ypix][xpix];
       if(pcr::isMV(val_x1y1))
       {
           tmv ++;
       }else {
           tvaltemp = val_x1y1;
        }
    }else {
       pcr::setMV(val_x1y1);
       tmv ++;
    }

    if(!((xpix_x < 0)|| (xpix_x > (nrCols() - 1)) || (ypix < 0)|| (ypix > (nrRows() - 1))))
    {
       val_x2y1 = data[ypix][xpix_x];
       if(pcr::isMV(val_x2y1))
       {
           tmv ++;
       }else {
           tvaltemp = val_x2y1;
        }
    }else {
       pcr::setMV(val_x2y1);
       tmv ++;
    }

    if(!((xpix < 0)|| (xpix > (nrCols() - 1)) || (ypix_y < 0)|| (ypix_y > (nrRows() - 1))))
    {
       val_x1y2 = data[ypix_y][xpix];
       if(pcr::isMV(val_x1y2))
       {
           bmv ++;
       }else {
           bvaltemp = val_x1y2;
        }
    }else {
       pcr::setMV(val_x1y2);
       bmv ++;
    }

    if(!((xpix_x < 0)|| (xpix_x > (nrCols() - 1)) || (ypix_y < 0)|| (ypix_y > (nrRows() - 1))))
    {
       val_x2y2 = data[ypix_y][xpix_x];
       if(pcr::isMV(val_x2y2))
       {
           bmv ++;
       }else {
           bvaltemp = val_x2y2;
        }
    }else {
       pcr::setMV(val_x2y2);
       bmv++;
    }

    if((tmv == 2) && (bmv == 2))
    {
        float val;
        pcr::setMV(val);
        return val;

    }else if((tmv == 1) && (bmv == 2))
    {
        return tvaltemp;
    }else if((tmv == 2) && (bmv == 1))
    {
        return bvaltemp;
    }else if((tmv == 0) && (bmv == 2))
    {
        return x_frac * val_x2y1 + (1.0-x_frac)* val_x1y1;
    }else if((tmv == 2) && (bmv == 0))
    {
        return x_frac * val_x1y2 + (1.0-x_frac)* val_x2y2;
    }else if((tmv == 1) && (bmv == 0))
    {
        return y_frac * ( x_frac * val_x1y2 + (1.0-x_frac)* val_x2y2) + (1.0-y_frac)* tvaltemp;
    }else if((tmv == 1) && (bmv == 1))
    {
        return y_frac * bvaltemp + (1.0-y_frac)* tvaltemp;
    }else if((tmv == 0) && (bmv == 1))
    {
        return y_frac * bvaltemp + (1.0-y_frac)* ( x_frac * val_x1y1 + (1.0-x_frac)* val_x2y1) ;
    }else if((tmv == 0) && (bmv == 0))
    {
        return y_frac * ( x_frac * val_x1y2 + (1.0-x_frac)* val_x2y2) + (1.0-y_frac)* ( x_frac * val_x1y1 + (1.0-x_frac)* val_x2y1) ;
    }

    return data[ypix][xpix];
}
template<typename T>
inline T cTMap<T>::GetPixelCovered(int r, int c, float cover)
{
    T ret = cover;

    if(r > -1 && r< nrRows() && c > -1 && c < nrCols())
    {
        T val = data[r][c];
        if(!pcr::isMV(val))
        {
            ret = val;
        }
    }
    return ret;
}
template<typename T>
inline T cTMap<T>::SampleCubic(double x, double y)
{
    // calculate coordinates -> also need to offset by half a pixel to keep image from shifting down and left half a pixel

    x = (x - 0.5 * cellSizeX() - west())/cellSizeX();
    y = (y - 0.5 * cellSizeY() - north())/cellSizeY();

    int xint = int(x);
    T xfract = std::fabs((x - floor(x)));//(cellSizeY() > 0?  std::fabs((x - floor(x))) : 1.0 - std::fabs((x - floor(x))));

    int yint = int(y);
    T yfract = std::fabs((y - floor(y)));// (cellSizeY() > 0? std::fabs( (y - floor(y))) : 1.0 - std::fabs( (y - floor(y))));
    T value;

    if(!((xint < 0)|| (xint > (nrCols() - 1)) || (yint < 0)|| (yint> (nrRows() - 1))))
    {

        value = data[yint][xint];

    }else {
        T x;
        pcr::setMV(x);
        return x;
    }

    int signx = 1.0;//(cellSizeX() > 0? 1:-1);
    int signy = 1.0;//(cellSizeY() > 0? 1:-1);

    // 1st row
    float p00 = GetPixelCovered(yint-signy*1,xint-signx*1,value);
    float p10 = GetPixelCovered(yint-signy*1,xint+signx*0,value);//GetPixelClamped(image, xint + 0, yint - 1);
    float p20 = GetPixelCovered(yint-signy*1,xint+signx*1,value);//GetPixelClamped(image, xint + 1, yint - 1);
    float p30 = GetPixelCovered(yint-signy*1,xint+signx*2,value);//GetPixelClamped(image, xint + 2, yint - 1);

    // 2nd row
    float p01 = GetPixelCovered(yint+signy*0,xint-signx*1,value);//GetPixelClamped(image, xint - 1, yint + 0);
    float p11 = GetPixelCovered(yint+signy*0,xint+signx*0,value);//GetPixelClamped(image, xint + 0, yint + 0);
    float p21 = GetPixelCovered(yint+signy*0,xint+signx*1,value);//GetPixelClamped(image, xint + 1, yint + 0);
    float p31 = GetPixelCovered(yint+signy*0,xint+signx*2,value);//GetPixelClamped(image, xint + 2, yint + 0);

    // 3rd row
    auto p02 = GetPixelCovered(yint+signy*1,xint-signx*1,value);//GetPixelClamped(image, xint - 1, yint + 1);
    auto p12 = GetPixelCovered(yint+signy*1,xint+signx*0,value);//GetPixelClamped(image, xint + 0, yint + 1);
    auto p22 = GetPixelCovered(yint+signy*1,xint+signx*1,value);//GetPixelClamped(image, xint + 1, yint + 1);
    auto p32 = GetPixelCovered(yint+signy*1,xint+signx*2,value);//GetPixelClamped(image, xint + 2, yint + 1);

    // 4th row
    float p03 = GetPixelCovered(yint+signy*2,xint-signx*1,value);//GetPixelClamped(image, xint - 1, yint + 2);
    float p13 = GetPixelCovered(yint+signy*2,xint+signx*0,value);//GetPixelClamped(image, xint + 0, yint + 2);
    float p23 = GetPixelCovered(yint+signy*2,xint+signx*1,value);//GetPixelClamped(image, xint + 1, yint + 2);
    float p33 = GetPixelCovered(yint+signy*2,xint+signx*2,value);//GetPixelClamped(image, xint + 2, yint + 2);

    // interpolate bi-cubically!

    float col0 = LISEMMath::CubicHermite(p00, p10, p20, p30, xfract);
    float col1 = LISEMMath::CubicHermite(p01, p11, p21, p31, xfract);
    float col2 = LISEMMath::CubicHermite(p02, p12, p22, p32, xfract);
    float col3 = LISEMMath::CubicHermite(p03, p13, p23, p33, xfract);
    float valuen = LISEMMath::CubicHermite(col0, col1, col2, col3, yfract);

    return valuen;

}


template<typename T>
inline T cTMap<T>::Sample(double x, double y)
{
    float ans = 0.0;

    float x_pix = (x - west())/cellSizeX();
    float y_pix = (y - north())/cellSizeY();

    //get closest pixel

    int xpix = x_pix;
    int ypix = y_pix;

    if((xpix < 0)|| (xpix > (nrCols() - 1)))
    {
        pcr::setMV(ans);
        return ans;
    }
    if((ypix < 0)|| (ypix > (nrRows() - 1)))
    {
        pcr::setMV(ans);
        return ans;
    }

    return data[ypix][xpix];
}


template<typename T>
inline T cTMap<T>::DerX( int r, int c)
{
    cTMap<T> *m =this;
    if(m != nullptr)
    {

        float val = 0.0;
        float val1 = 0.0;
        float val2 = 0.0;
        if(r > 0 && c > 0 && r < m->nrRows() && c < m->nrCols())
        {
            val = m->data[r][c];
            if(pcr::isMV(val))
            {
                val = 0.0;
            }
        }
        if(r > 0 && c-1 > 0 && r < m->nrRows() && c-1 < m->nrCols())
        {
            val1 = m->data[r][c-1];
            if(pcr::isMV(val1))
            {
                val1 = val;
            }
        }
        if(r > 0 && c+1 > 0 && r < m->nrRows() && c+1 < m->nrCols())
        {
            val2 = m->data[r][c+1];
            if(pcr::isMV(val2))
            {
                val2 = val;
            }
        }
        return 0.5 * ((val2 - val) + (val - val1));
    }
    return 0.0;
}
template<typename T>
inline T cTMap<T>::DerY( int r, int c)
{
    cTMap<T> *m =this;
    if(m != nullptr)
    {

        float val = 0.0;
        float val1 = 0.0;
        float val2 = 0.0;
        if(r > 0 && c > 0 && r < m->nrRows() && c < m->nrCols())
        {
            val = m->data[r][c];
            if(pcr::isMV(val))
            {
                val = 0.0;
            }
        }
        if(r-1 > 0 && c > 0 && r-1 < m->nrRows() && c < m->nrCols())
        {
            val1 = m->data[r-1][c];
            if(pcr::isMV(val1))
            {
                val1 = val;
            }
        }
        if(r+1 > 0 && c > 0 && r+1 < m->nrRows() && c < m->nrCols())
        {
            val2 = m->data[r+1][c];
            if(pcr::isMV(val2))
            {
                val2 = val;
            }
        }
        return 0.5 * ((val2 - val) + (val - val1));
    }
    return 0.0;
}

template<typename T>
inline bool cTMap<T>::CheckIsEdgeOfMaskedDomain(size_t r, size_t c){
    if(r == 0 || c == 0 || r  == nrRows()- 1 || c == nrCols() -1)
    {
        if(!pcr::isMV(data[r][c]))
        {
            return true;
        }else
        {
            return false;
        }
    }
    if(pcr::isMV(data[r][c]))
    {
        return false;
    }

    if(pcr::isMV(data[std::max((size_t)0,std::min((size_t)nrRows(),r))][std::max((size_t)0,std::min((size_t)nrCols(),c+1))]))
    {
        return true;
    }
    if(pcr::isMV(data[std::max((size_t)0,std::min((size_t)nrRows(),r))][std::max((size_t)0,std::min((size_t)nrCols(),c-1))]))
    {
        return true;
    }
    if(pcr::isMV(data[std::max((size_t)0,std::min((size_t)nrRows(),r+1))][std::max((size_t)0,std::min((size_t)nrCols(),c))]))
    {
        return true;
    }
    if(pcr::isMV(data[std::max((size_t)0,std::min((size_t)nrRows(),r-1))][std::max((size_t)0,std::min((size_t)nrCols(),c))]))
    {
        return true;
    }
    return false;
}
template<typename T>
inline std::string cTMap<T>::projectionstring() const
{
    return m_projection;
}
template<typename T>
inline void cTMap<T>::setprojectionstring(std::string proj)
{
    m_projection = proj;
}
template<typename T>
inline BoundingBox cTMap<T>::GetExtent()
{
    return BoundingBox(west(),west() +cellSizeX()*((double)(nrCols())),north(), north() +cellSizeY()*((double)(nrRows())));
}




template<typename T>
inline void DeNormalizeMapInPlace(cTMap<T> *m, double min, double max)
{

    bool found = false;
    //get minimum and maximum
    for(int r = 0; r < m->data.nr_rows();r++)
    {
        for(int c = 0; c < m->data.nr_cols();c++)
        {
            if(!pcr::isMV(m->data[r][c]))
            {
                m->data[r][c] = min + (max-min) * m->data[r][c];
            }
        }
    }
}
template<typename T>
inline void NormalizeMapInPlace(cTMap<T> *m, double &min, double &max)
{

    bool found = false;
    //get minimum and maximum
    for(int r = 0; r < m->data.nr_rows();r++)
    {
        for(int c = 0; c < m->data.nr_cols();c++)
        {
            if(!pcr::isMV(m->data[r][c]))
            {
                if(!found)
                {
                    min = m->data[r][c];
                    max = m->data[r][c];
                    found = true;
                }else {

                    min = std::min(min,(double)m->data[r][c]);
                    max = std::max(max,(double)m->data[r][c]);
                }
            }
        }
    }
    //normalize values between zero and 1
    if(max > min)
    {
        for(int r = 0; r < m->data.nr_rows();r++)
        {
            for(int c = 0; c < m->data.nr_cols();c++)
            {
                if(!pcr::isMV(m->data[r][c]))
                {
                    m->data[r][c] = std::min(1.0,std::max(0.0,(m->data[r][c] - min)/(max-min)));
                }else {
                    pcr::setMV(m->data[r][c]);
                }
            }
        }
    }
}

template<typename T>
inline cTMap<T> * NormalizeMap(cTMap<T> *m, double &min, double &max)
{
    cTMap<T> * res = m->GetCopy();

    bool found = false;
    //get minimum and maximum
    for(int r = 0; r < m->data.nr_rows();r++)
    {
        for(int c = 0; c < m->data.nr_cols();c++)
        {
            if(!pcr::isMV(m->data[r][c]))
            {
                if(!found)
                {
                    min = m->data[r][c];
                    max = m->data[r][c];
                    found = true;
                }else {

                    min = std::min(min,(double)m->data[r][c]);
                    max = std::max(max,(double)m->data[r][c]);
                }
            }
        }
    }
    //normalize values between zero and 1
    if(max > min)
    {
        for(int r = 0; r < m->data.nr_rows();r++)
        {
            for(int c = 0; c < m->data.nr_cols();c++)
            {
                if(!pcr::isMV(m->data[r][c]))
                {
                    res->data[r][c] = std::min(1.0,std::max(0.0,(res->data[r][c] - min)/(max-min)));
                }
            }
        }
    }

    return res;
}

template<typename T>
inline bool OUTORSMALLER(cTMap<T> * map, int r, int c, double crit)
{
    if(r<0 || c< 0|| r>map->nrRows()-1 || c > map->nrCols()-1)
    {
        return true;
    }
    if(!pcr::isMV(map->data[r][c]) && map->data[r][c] < crit)
    {
        return true;
    }
    return false;
}


template<typename T>
inline bool OUTORMV(cTMap<T> * map, int r, int c)
{
    if(r<0 || c< 0|| r>map->nrRows()-1 || c > map->nrCols()-1)
    {
        return true;
    }
    if(pcr::isMV(map->data[r][c]))
    {
        return true;
    }
    return false;
}

template<typename T>
inline bool OUTORMV(std::shared_ptr<cTMap<T>> map, int r, int c)
{
    if(r<0 || c< 0|| r>map->nrRows()-1 || c > map->nrCols()-1)
    {
        return true;
    }
    if(pcr::isMV(map->data[r][c]))
    {
        return true;
    }
    return false;
}

template<typename T>
LSMVector2  cTMap<T>::GetPixelLngLat       (int r, int c)
{
    GeoProjection pfrom = GeoProjection::FromString(this->projectionstring());
    GeoProjection pto = GeoProjection::GetFromEPSG(4326);
    std::shared_ptr<GeoCoordTransformer>tr = GeoCoordTransformer::GetCoordTransformerSP(pfrom,pto);

    if(tr)
    {
        return tr->Transform(LSMVector2(this->west() + c * this->cellSizeX(),this->north() + r * this->cellSizeY()));
    }else
    {
        return LSMVector2(this->west() + c * this->cellSizeX(),this->north() + r * this->cellSizeY());
    }

}

//===========================================================
//          Map Arithmatic
//===========================================================


//forward declaration of map arithmatic base function
template<typename T, typename Callable>
cTMap<T> MapArithmatic(const cTMap<T> &map1,const cTMap<T> &map2,const Callable &f)
{
    if(!(map1.nrRows() == map2.nrRows() && map1.nrCols() == map2.nrCols()))
    {
        throw std::runtime_error("Error for raster arithmatic, mismatch rows and columns; m1: " + std::to_string(map1.nrCols()) +  "x" + std::to_string(map1.nrRows()) +" m2: "+ std::to_string(map2.nrCols()) +  "x" + std::to_string(map2.nrRows()));
    }
    cTMap<T> ret = map1;

#if(FF_OMP)
    #pragma omp parallel for collapse(2)
    for(int r = 0; r < map1.data.nr_rows();r++)
    {
        for(int c = 0; c < map1.data.nr_cols();c++)
        {
            if(pcr::isMV(map1.data[r][c]) || pcr::isMV(map2.data[r][c]) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(map1.data[r][c],map2.data[r][c]);
            }
        }
    }

#else
    GetTPool().parallel_for(0,map1.data.nr_rows(),[&](int r)
    {
        for(int c = 0; c < map1.data.nr_cols();c++)
        {
            if(pcr::isMV(map1.data[r][c]) || pcr::isMV(map2.data[r][c]) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(map1.data[r][c],map2.data[r][c]);
            }
        }

    });

    GetTPool().waitext();

#endif

    return ret;
}


//forward declaration of map arithmatic base function
template<typename T, typename G, typename Callable>
cTMap<T> MapRealArithmatic(const cTMap<T> &map1,const G &map2,const Callable &f)
{
    cTMap<T> ret = map1;

#if(FF_OMP)
    #pragma omp parallel for collapse(2)
    for(int r = 0; r < map1.data.nr_rows();r++)
    {
        for(int c = 0; c < map1.data.nr_cols();c++)
        {
            if(pcr::isMV(map1.data[r][c]) || pcr::isMV(map2) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(map1.data[r][c],map2);
            }
        }
    }

#else
    GetTPool().parallel_for(0,map1.data.nr_rows(),[&](int r)
    {
        for(int c = 0; c < map1.data.nr_cols();c++)
        {
            if(pcr::isMV(map1.data[r][c]) || pcr::isMV(map2) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(map1.data[r][c],map2);
            }
        }

    });
    GetTPool().waitext();

#endif

    return ret;
}

//forward declaration of map arithmatic base function
template<typename T, typename Callable>
void MapArithmaticInPlace(cTMap<T> &map1,const cTMap<T> &map2,const Callable &f)
{
    if(!(map1.nrRows() == map2.nrRows() && map1.nrCols() == map2.nrCols()))
    {
        throw std::runtime_error("Error for raster arithmatic, mismatch rows and columns; m1: " + std::to_string(map1.nrCols()) +  "x" + std::to_string(map1.nrRows()) +" m2: "+ std::to_string(map2.nrCols()) +  "x" + std::to_string(map2.nrRows()));
    }
    cTMap<T> &ret = map1;

#if(FF_OMP)
    #pragma omp parallel for collapse(2)
    for(int r = 0; r < map1.data.nr_rows();r++)
    {
        for(int c = 0; c < map1.data.nr_cols();c++)
        {
            if(pcr::isMV(map1.data[r][c]) || pcr::isMV(map2.data[r][c]) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(map1.data[r][c],map2.data[r][c]);
            }
        }
    }

#else
    GetTPool().parallel_for(0,map1.data.nr_rows(),[&](int r)
    {
        for(int c = 0; c < map1.data.nr_cols();c++)
        {
            if(pcr::isMV(map1.data[r][c]) || pcr::isMV(map2.data[r][c]) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(map1.data[r][c],map2.data[r][c]);
            }
        }

    });
    GetTPool().waitext();

#endif

}


//forward declaration of map arithmatic base function
template<typename T, typename Callable>
void MapRealArithmaticInPlace(cTMap<T> &map1,const T &map2,const Callable &f)
{
    cTMap<T> &ret = map1;

#if(FF_OMP)
    #pragma omp parallel for collapse(2)
    for(int r = 0; r < map1.data.nr_rows();r++)
    {
        for(int c = 0; c < map1.data.nr_cols();c++)
        {
            if(pcr::isMV(map1.data[r][c]) || pcr::isMV(map2) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(map1.data[r][c],map2);
            }
        }
    }

#else
    GetTPool().parallel_for(0,map1.data.nr_rows(),[&](int r)
    {
        for(int c = 0; c < map1.data.nr_cols();c++)
        {
            if(pcr::isMV(map1.data[r][c]) || pcr::isMV(map2) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(map1.data[r][c],map2);
            }
        }

    });
    GetTPool().waitext();

#endif

}

template<typename T,  typename Callable>
cTMap<T> CallMapSingleFunc(const cTMap<T> & m,const Callable &f)
{
    cTMap<T> ret = m;

#if(FF_OMP)
    #pragma omp parallel for collapse(2)
    for(int r = 0; r < m.data.nr_rows();r++)
    {
        for(int c = 0; c < m.data.nr_cols();c++)
        {
            if(pcr::isMV(m.data[r][c]) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(m.data[r][c]);
            }
        }
    }

#else
    GetTPool().parallel_for(0,m.data.nr_rows(),[&](int r)
    {
        for(int c = 0; c < m.data.nr_cols();c++)
        {
            if(pcr::isMV(m.data[r][c]) )
            {
                pcr::setMV(ret.data[r][c]);
            }else
            {
                ret.data[r][c] = f(m.data[r][c]);
            }
        }

    });
    GetTPool().waitext();

#endif

    return ret;
}

template <typename T>
inline cTMap<T>* MakeMapCopy(const cTMap<T> * m, T val)
{
    cTMap<T> * m2= m->GetCopy();
    m2->setAllValue(val);
    return m2;
}



inline void Safe_Delete_MapVector(std::vector<cTMapf*>v)
{
    for(int i = 0; i < v.size(); i++)
    {
        if(!(v[i] == nullptr))
        {
            delete v[i];
        }
    }


}

#endif
