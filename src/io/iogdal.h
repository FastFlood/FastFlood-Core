#ifndef IOGDAL_H
#define IOGDAL_H

#include "map.h"
#include <string>
#include <iostream>
#include <filesystem>

#if FF_HAS_GDAL
#include "gdal.h"
#include "gdal_priv.h"
#include "cpl_string.h"
#include "cpl_conv.h"
#include "cpl_vsi.h"
#include "ogr_api.h"
#include "ogrsf_frmts.h"
#include "gdalwarper.h"
#include "gdal_alg.h"
#include "gdal_utils.h"
#endif


#if FF_HAS_GDAL
inline static bool PreInitializeGDAL(std::string working_path)
{
    return true;
}

//! Function to close a GDAL GDALDataset.
const auto close_gdal_dataset = [](GDALDataset* dataset) { GDALClose(dataset); };

//! Auto-ptr type for GDAL GDALDatasets.
using GDALDatasetPtr = std::unique_ptr<GDALDataset, decltype(
    close_gdal_dataset)>;

inline static void GDALERROR(CPLErr er, CPLErrorNum ernum, const char* message)
{
    std::cout << "Error from gdal (" << ernum << ") " << message << std::endl;
}


inline static bool InitializeGDAL()
{
    CPLSetConfigOption("GDAL_MAX_BAND_COUNT", "1000000");
    CPLSetConfigOption("HDF5_DISABLE_VERSION_CHECK", "1");

    CPLSetErrorHandler(GDALERROR);

    GDALAllRegister();
    OGRRegisterAll();

    return true;
}

//Read the cTMapf raster structure using GDAL
inline cTMapf *readRaster(
    std::string const& pathName, int bandn = 0)
{

    std::cout << "Read raster " <<  pathName << "  Bandnr: " << bandn << std::endl;

    // Open raster dataset and obtain some properties.
    GDALDataset* dataset = static_cast<GDALDataset*>(GDALOpen(
        pathName.data(), GA_ReadOnly));
    if(!dataset) {
        //ff_error = "Can not open file at location: " + pathName;
        return nullptr;
    }

    int nr_bands = dataset->GetRasterCount();
    if(nr_bands == 0) {
        //ff_error = "Raster file has no bands: " + pathName;
        return nullptr;
    }

    double transformation[6];
    dataset->GetGeoTransform(transformation);
    const OGRSpatialReference * ref = dataset->GetSpatialRef();
    char * refc;
    if(ref != nullptr)
    {
        ref->exportToWkt(&refc);
    }

    std::string projection{dataset->GetProjectionRef()};


    // Read the first raster band.
    GDALRasterBand* band{dataset->GetRasterBand(1 + bandn)};
    assert(band);

    int const nr_rows{band->GetYSize()};
    int const nr_cols{band->GetXSize()};
    double const west{transformation[0]};
    double const north{transformation[3]};
    double const cell_size{transformation[1]};
    double const cell_sizeY{transformation[5]};


    MaskedRaster<float> raster_data(nr_rows, nr_cols, north, west, cell_size, cell_sizeY);


    if(band->RasterIO(GF_Read, 0, 0, nr_cols, nr_rows, raster_data[0],
            nr_cols, nr_rows, GDT_Float32, 0, 0) != CE_None) {
        //ff_error = "Raster file has no bands: " + pathName;
        return nullptr;
    }

    int hasNoDataValue{false};
    double noDataValue{band->GetNoDataValue(&hasNoDataValue)};
    if(hasNoDataValue) {
        raster_data.replace_with_mv(noDataValue);
    }

    GDALClose((GDALDatasetH) dataset);

    cTMapf ret1 = cTMapf(std::move(raster_data), projection,pathName,false);

    cTMapf * ret = new cTMapf(std::move(ret1));
    return ret;
}


inline int writeGDALRaster(
    cTMapf const& raster,
    std::string const& pathName,
    GDALDriver& driver,
    int bandn = -1, bool compress = false)
{

    CPLStringList sl;
    if(compress)
    {
        sl.AddNameValue("COMPRESS","DEFLATE");
        sl.AddNameValue("PREDICTOR","3");
        sl.AddNameValue("Tiled","YES");
        sl.AddNameValue("BLOCKXSIZE","256");
        sl.AddNameValue("BLOCKYSIZE","256");
    }
    try
    {



    if(bandn == -1)
    {
        if(std::filesystem::exists(pathName))
        {
            std::cout << "Output file already exists (try and remove): " + pathName << std::endl;
            int remove_code = remove(pathName.data());
            if(remove_code != 0)
            {
                std::cout << "Could not remove file" <<std::endl;
            }
        }
        int const nrRows{raster.nrRows()};
        int const nrCols{raster.nrCols()};
        int const nrBands{1};
        GDALDataset * d;

        d = driver.Create(pathName.data(),
                      nrCols, nrRows, nrBands, GDT_Float32, sl.List());

        if(d == NULL) {
             //ff_error = "Dataset cannot be created: " + pathName;

             const char * c = CPLGetLastErrorMsg();
             if(c!=nullptr)
             {
                 //ff_error = "Dataset cannot be created: " + pathName + " error: " + std::string(c);
             }
             return -1;
        }


        MaskedRaster<float> const& raster_data{raster.data};

        double transformation[]{
            raster_data.west(),
            raster_data.cell_size(),
            0.0,
            raster_data.north(),
            0.0,
            raster_data.cell_sizeY()};

        d->SetGeoTransform(transformation);


        if(!(raster.projectionstring() == "LISEM_GENERIC" || raster.projectionstring() == ""))
        {
            d->SetProjection(raster.projectionstring().data());
        }

        d->SetMetadataItem("PCRASTER_VALUESCALE", "VS_SCALAR");

        auto band = d->GetRasterBand(1);

        band->SetNoDataValue(-FLT_MAX);

        if(band->RasterIO(GF_Write, 0, 0, nrCols, nrRows,
                const_cast<float*>(&raster_data.cell(0)),
                nrCols, nrRows, GDT_Float32, 0, 0) != CE_None) {
            //ff_error = "Raster band cannot be written: " + pathName;
            return -1;
        }


        GDALClose( (GDALDataset*) d );


        return 0;
    }else {

        std::cout << "try to append to existing file";

        GDALDataset* dataset(static_cast<GDALDataset*>(GDALOpen(
            pathName.data(), GA_Update)));
        if(!dataset) {

            int const nrRows{raster.nrRows()};
            int const nrCols{raster.nrCols()};
            int const nrBands{1};
            GDALDataset * d = driver.Create(pathName.data(),
                          nrCols, nrRows, nrBands, GDT_Float32, nullptr);

            if(d == NULL) {
                 //ff_error = "Dataset cannot be created: " + pathName;
                 return -1;
            }


            MaskedRaster<float> const& raster_data{raster.data};

            double transformation[]{
                raster_data.west(),
                raster_data.cell_size(),
                0.0,
                raster_data.north(),
                0.0,
                raster_data.cell_sizeY()};

            d->SetGeoTransform(transformation);



            if(!(raster.projectionstring() == "LISEM_GENERIC" || raster.projectionstring() == ""))
            {
                d->SetProjection(raster.projectionstring().data());
            }

            d->SetMetadataItem("PCRASTER_VALUESCALE", "VS_SCALAR");

            auto band = d->GetRasterBand(1);

            band->SetNoDataValue(-FLT_MAX);


            std::cout << "writing to band create"<< band->GetAccess() << std::endl;

            if(band->RasterIO(GF_Write, 0, 0, nrCols, nrRows,
                    const_cast<float*>(&raster_data.cell(0)),
                    nrCols, nrRows, GDT_Float32, 0, 0) != CE_None) {

                //ff_error = "Band can not be created: " + pathName;
                return -1;
            }


            GDALClose( (GDALDataset*) d );


            return 0;
        }

        int nr_bands = dataset->GetRasterCount();
        if(nr_bands < bandn) {

            //ff_error = "Map %1 does not contain band: " + pathName + " band nr: ";
            return -1;
        }

        auto band = dataset->GetRasterBand(1+bandn);

        band->SetNoDataValue(-FLT_MAX);

        int const nrRows{raster.nrRows()};
        int const nrCols{raster.nrCols()};

        MaskedRaster<float> const& raster_data{raster.data};

        if(band->RasterIO(GF_Write, 0, 0, nrCols, nrRows,
                const_cast<float*>(&raster_data.cell(0)),
                nrCols, nrRows, GDT_Float32, 0, 0) != CE_None) {
            //ff_error = "Band can not be written to file: " + pathName;
            return -1;

        }

        GDALClose( (GDALDataset*) dataset );

        return 0;


   }
    } catch (const std::exception& ex) {
        std::cout << "Error during writing of georaster: " << ex.what() <<std::endl;
    } catch (const std::string& ex) {
        std::cout << "Error during writing of georaster: " << ex <<std::endl;
    }catch(...)
    {
        std::cout << "Error during writing of georaster: " << "unkown" <<std::endl;
    }
    return 0;
}

inline int writeRaster(
    cTMapf const& raster,
    std::string const& pathName,
    std::string format = "GTIFF",
    int bandn = -1,
        bool compress = true)
{
    if(raster.nrRows() == 0 || raster.nrCols() == 0) {
        return -1;
    }

    if(pathName == "") {

        return -1;
    }

    GDALDriver* driver = GetGDALDriverManager()->GetDriverByName(
        format.data());

    if(!driver) {

        //ff_error = "Raster driver does not exist " + format;
        return -1;
    }

    char** metadata{driver->GetMetadata()};
    bool driverSupportsCreate{CSLFetchBoolean(metadata, GDAL_DCAP_CREATE,
        FALSE) != FALSE};

    if(driverSupportsCreate) {
        return writeGDALRaster(raster, pathName, *driver,bandn,compress);
    }else {
        //ff_error = "Raster driver can not be used to create dataset " + format;
        return -1;
    }

    return 0;
}


inline cTMapf *load_map(std::string file)
{
    cTMapf * m = readRaster(file,0);

    return m;
}

#endif
#endif // IOGDAL_H
