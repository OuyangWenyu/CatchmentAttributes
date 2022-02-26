# CatchmentAttributes

This is a Catchment-Attributes Calculator (Mainly for Contiguous China now).

This repository is forked from https://github.com/haozhen315/CCAM-China-Catchment-Attributes-and-Meteorology-dataset

The original repository contains the complement code for paper:

```
Hao, Z., Jin, J., Xia, R., Tian, S., Yang, W., Liu, Q., Zhu, M., Ma, T., Jing, C., and Zhang, Y.: CCAM: China Catchment Attributes and Meteorology dataset, Earth Syst. Sci. Data, 13, 5591–5616, https://doi.org/10.5194/essd-13-5591-2021, 2021.
```

The manuscript can be found [here](https://doi.org/10.5194/essd-13-5591-2021) (publicly available).

Now [I, Wenyu Ouyang](https://github.com/OuyangWenyu), modify it to an attribute calculator for catchments. It supports
generating 120+ basin attributes for each basin given a single or several basin boundaries. Except for the netCDF data,
which only covers China, other data sources are global. Now all catchments should be divided up to the source of the
river such that there is no upper catchment.

In the following sections, I will use CatchmentAttributes to represent my version, and original CatchmentAttributes or
original repository for the previous version.

If you wanna use CatchmentAttributes to calculate your basins' attributes, please read this README.md carefully and use
it step by step. In the future, we may provide a more intelligent version.

If you find any problems or unclear in the code, you can contact Wenyu Ouyang through wenyuouyang@outlook.com

## Prerequisite

### Shapefiles

This is a basins' attributes calculator. Hence, basins are the prerequisite.

Generally we have shapefiles for our target basins, and the outlets of basins are always streamflow stations.

But if you don't have shapefiles, you have to delineate basins by yourself before running CatchmentAttributes.

If you have little experience on delineating basins and can read Chinese well, please
see [this instruction](https://github.com/OuyangWenyu/hydroGIS/blob/dev/QGIS/2-watershed-delineation.md) for watershed
delineation.

"One shapefile for one basin" is required here. In addition, all shapefiles of basins should have a same coordination
system (we chose **EPSG:4326 - WGS 84**). To be convenient, please set the name of gage_id field as **"ID"**. After
getting the shapefiles, please put it in the "data/shapefiles" directory in the root directory of CatchmentAttributes.

Please give your basins ids so that we can process them in a uniform way. The template is "basin_xxx" (xxx means some
numbers). See an example in "data/shapefiles" directory.

Don't put any shapefiles in the "data/shapefiles" directory except for basins' shapefiles.

### QGIS

QGIS is used in many parts of CatchmentAttributes, so you'd better know some basic knowledge about it.
[Here](https://github.com/OuyangWenyu/hydroGIS/tree/dev/QGIS) is a very simple Chinese tutorial.

You can create a QGIS project in the "data" directory of CatchmentAttributes and use it to handle with all processes
during using CatchmentAttributes.

### NASA Earthdata

Some source data should be downloaded from NASA Earthdata. If you are new to it and can read Chinese, please
see [this instruction](https://github.com/OuyangWenyu/aqualord/tree/master/DEM) and try to download some data.

### Google Earth Engine

We will use [Google earth engine, GEE](https://code.earthengine.google.com/) to process the Meteorological forcing data.

You can learn GEE from [here](https://github.com/OuyangWenyu/hydroGIS/tree/dev/GEE) if you never used it before.

### Python

Of course, this is a python repository, so you have to know how to use python. We also provide a
Chinese [tutorial](https://github.com/waterDLut/hydrus/tree/dev) for modeling with python.

Quick steps to create python environment for CatchmentAttributes are as follows:

```Shell
# After downloading CatchmentAttributes from github, move to its root directory
cd CatchmentAttributes
# create conda environment
conda env create -f environment.yml
# activate the env
conda activate CatchAttr
```

## Meteorological time series

It seems that the situ observations meteorological data are no longer available in the data website of China
Meteorological Agency.

Even though you can find the source data, their downloading is very difficult especially you are a student rather than a
teacher.

Hence, we recommend to use [ERA5-LAND reanalysis data](https://essd.copernicus.org/articles/13/4349/2021/) to run
CatchmentAttributes. But notice the resolution of ERA5-LAND is about 10km. If your basins are very small, such as <
100km2, maybe this is not good choice.

1. Upload the shapefiles of your target basins (put all basins' polygons into one shapefile) to GEE and
   use [this script](https://code.earthengine.google.com/b40bfe7529ec9df928f5abfc029f4e4e) to calculate basin mean
   forcing values for each basin.
2. Then download the forcing data files from your Google Drive and unzip it. Change the unzipped directory's name to "
   ERA5_LAND"
3. run the following script to transform them to a standard format to make them easy-to-use.

```Shell
cd catch_attr
# Users should modify year_range and region to their own. 
# For example, mine data's time range is 2014-01-01 to 2020-12-31 
# and my source data file's name is era5_land_camels_mr_avg_mean_2010.csv where "camels_mr" is the region
python app.py --catch_attr basin_mean_forcing --year_range 2014 2021 --region camels_mr
```

Then, the transformed data will be in "data/basin_mean_forcing".

The original CatchmentAttributes repository provide methods for interpolating site observation climate data to rasters (
GeoTIFF). But we don't need them now. If you want to use its tools for other use, please
see [the original tutorial](https://github.com/haozhen315/catchment-attributes-and-meteorology-for-large-sample-study-in-contiguous-china)
and the code are still reserved here in catch_attr/raster_surf.py and catch_attr/raster2catchment.py (but we didn't test
them in CatchmentAttributes)

## Climate indicator

After getting forcing data, we can calculate the climate indicators.

Run the following code:

```Shell
python app.py --catch_attr climate
```

The result will be stored in data/attribute directory -- climate.csv

## Lithology

1. Download the GLiM dataset from original
   authors: https://www.dropbox.com/s/9vuowtebp9f1iud/LiMW_GIS%202015.gdb.zip?dl=0
2. Create a "GLIM" directory in the "data" directory, put the downloaded zip file in it and unzip the file to LiMW_GIS
   2015.gdb
3. Import the gdb dataset to QGis (use the QGIS project we build in Section "Prerequisite"): "Layer" -> "Add Layer" -> "
   Add Vector Layer", then choose the directory "LiMW_GIS 2015.gdb". The aim is to convert GliM to a raster form for
   processing.
4. If we directly deal with the whole file, it will cost too much; now that we have shapefiles, it's a better choice to
   clip the whole file to a range. You can use "Vector overlay" -> "Extract/clip by extent" and chose an extent which
   overlay all target basins to clip. If it didn't work, you can use "Vector geometry" -> "Fix geometries" to fix it at
   first. Remember all processed files are saved in the "GLIM" directory.
5. Export GLiM to GeoTIFF format. First, we need to use "add unique value index field" for the "xx" field ("xx" is the
   geology type) of the glim fixed shapefile to get new "xxValue" number field and output "layer with index field"
   to "xx_index_field.shp"; then, we can rasterize "xx_index_field.shp" ("xxValue") to a geotiff file and name it
   "glim_clip_raster.tif"
6. Run the following script:

```Shell
python app.py --catch_attr geology
```

The resulting file will appear in the data/attribute directory -- geology.csv

The files needed in the process:

```bash
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ......
├── GLiM
|   ├── LiMW_GIS 2015.gdb
|   ├── glim_clip.shp
|   ├── glim_clip_raster.tif
|   ├── glim_fix.shp
|   ├── GLiMCateNumberMapping.csv
|   ├── GlimRaster.tif
|   ├── LiMW_GIS 2015.gdb.zip
|   ├── xx_index_field.shp
├── glim_short_long_name.txt
```

The file "glim_name_short_long.txt" is not used in the process, but it is useful for understanding the geology attribute
file.

The file "data/glim_cate_number_mapping.csv" comes from original CatchmentAttributes, we didn't use it.

## Permeability and porosity

The [9_code_data.zip](https://zenodo.org/record/5137288/files/9_code_data.zip?download=1) included in the Zenodo
repository contains processed_glim.py, processed_igbp.tif, processed_permeability.tif and processed_porosity.tif

Here we need processed_permeability.tif and processed_porosity.tif.

Put the downloaded zip file in the "data" directory and unzip it.

The files needed in the process:

```bash
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ......
├── 9_code_data
|   ├── processed_permeability.tif
|   ├── processed_porosity.tif
|   ├── ...
```

Then you can run the following code:

```Shell
python app.py --catch_attr permeability_porosity
```

## Land cover

1. Source data: https://lpdaac.usgs.gov/products/mcd12q1v006/. However, MODIS data is divided into different tiles,
   which is inconvenient for processing. Original authors have merged the MODIS product into a single tif which can be
   downloaded here: https://1drv.ms/u/s!AqzR0fLyn9KKspF4xxbe0xM7qJNzkA?e=TYyZeC. Download it (processed_igbp.tif) and
   put it in the "data" directory
2. Run the following script:

```Shell
python app.py --catch_attr land_cover
```

The resulting file will appear in the data/attribute directory -- land_cover.csv.

Required data folder structure:

```bash
(1) processed_igbp.tif: converted IGBP classification in raster form
(2) Catchment shapefiles
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ...
```

## Root depth

Run the following script:

```Shell
python app.py --catch_attr root_depth
```

Required folder structure:

```bash
(1) processed_igbp.tif: converted IGBP classification in raster form
(2) calculated_root_depth.txt: calculated root_depth 50/99 for each type of land cover based on Eq. (2) and Table 2 in (Zeng 2001)
(2) Catchment shapefiles
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ...
```

## Topography (Elev and Slope)

1. Download ASTER GDEM from NASA Earthdata and put them in the "data/DEM" directory.
2. Run the following script:

```Shell
python app.py --catch_attr topo_elev
```

This will take long time: about a few minutes for one basin.

Required folder structure:

```bash
(1) ASTER GDEM
├── DEM
|   ├── ASTGTMV003_NxxExxx_dem.tif
|   ├── ...
(2) Catchment shapefiles
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ...
```

## Topography (Shape characteristics)

Shape characteristics are calculated based on catchment length (the mainstream length measured from the basin outlet to
the remotest point on the basin boundary), catchment perimeters and catchment area. The calculation is first performed
in decimal degree, and then the units are converted to kilometres by projecting to UTM coordinate system.

1. Download Asian river network data from: https://www.cger.nies.go.jp/db/gdbd/gdbd_index_e.html. Unzip it and import
   the .mdb file to QGIS. If you don't know how to import .mdb file, please refer
   to [this tutorial](https://github.com/OuyangWenyu/hydroGIS/blob/master/QGIS/3-qgis-common.md#%E5%AF%BC%E5%85%A5%E6%95%B0%E6%8D%AE)
2. The needed file is "as_streams". Now we can export from the .mdb database to a shp file ("data/TOPO/as_streams.shp")
   using QGIS (project it to WGS84 when exporting).
3. Run the following script

```Shell
python app.py --catch_attr topo_shape
```

Because as_streams.shp only provided some large rivers, for some basins, the code cannot find the outlet of the basins.
Hence, some basins' shape characteristics are NaN values.

Some warnings will occur, but they have no impact on the results; we will fix them in the future.

## LAI/NDVI

1. Download modis product: https://lpdaac.usgs.gov/products/mcd15a3hv006/ (for LAI, search "mcd15a3 006" in Earthdata
   website) and https://lpdaac.usgs.gov/products/mod13q1v006/ (for NDVI, search "mod13q1 006" in Earthdata website) from
   NASA Earthdata. The original CatchmentAttributes use 20 years' data from 2000-2020. If you think it's too long to
   download, you can chose a few years. For example, I only use data in 2015 for my target region (the files are still
   dozens of GBs).
2. The source data are in hdfs format. The provided script first find needed hdfs tiles for the given catchment and
   merge them. Then perform zonal statistics to get catchment-averaged values. Put the downloaded hdfs files into the
   folder "data/MCD15A3[MOD13Q1]" and create an output folder e.g. "data/LAI[NDVI]", and run the following script:

```Shell
# for LAI
python app.py --catch_attr modis --modis_type lai --year_range 2015 2016
# for NDVI
python app.py --catch_attr modis --modis_type ndvi --year_range 2015 2016
```

Required folder structure:

```bash
(1) MODIS data
├── MOD13Q1/MCD15A3H
|   ├── MCD15A3H.A2002185.h22v04.006.2015149102803.hdf
|   ├── MCD15A3H.A2002186.h22v04.006.2015149102803.hdf
|   ├── MCD15A3H.A2002187.h22v04.006.2015149102803.hdf
|   ├── ......
(2) Catchment shapefiles
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ...
```

## Soil

1. Download source files from [SoilGrids250](https://files.isric.org/soilgrids/former/2017-03-10/data/). As the files
   are large, you can only choose the attributes you want. More details about the attributes could be
   seen [here](https://github.com/ISRICWorldSoil/SoilGrids250m/blob/master/grids/models/META_GEOTIFF_1B.csv)
2. Downscale the downloaded tif files, so that we can process it later easily. Open a tif file in QGIS, then choose
   Raster->Projections->Warp (Reproject), and name the output file as "xxx_downscaled.tif". Notice the unit of
   resolution in this setting is same with the original tif file's and it is degree as the CRS is WGS84 (here we choose
   0.05). Because we will calculate the mean value of grids in a basin, here we also choose "Average" as the resampling
   method. Set no data value as -9999
3. put the processed tif files to "data/soil_source_data" directory. Then you can run the following code

```Shell
python app.py --catch_attr soil
```