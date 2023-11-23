import os
import re
import numpy as np
from osgeo import gdal, osr

import fiona
import rasterio
import rasterio.mask
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling


def geotif_from_array(
    array: np.array,
    lat_start: float,
    lat_end: float,
    lon_start: float,
    lon_end: float,
    degree: float,
    output_file: str,
):
    """
    将一个 numpy array 写入一个带有位置信息的 tif 文件，默认使用 wgs84 坐标系

    >>> geotif_from_array(array=res, lat_start=19.94174, lat_end=49.18826, lon_start=75.46174, lon_end=130.5756, degree=0.11652 ,output_file='results/tmp.tif')

    Parameters
    ----------
    array
        要写入 tif 的变量，其 shape 应该和 lats 和 lons 对应
    lat_start
        起始维度，可参考 arcmap 生成的栅格文件 source 属性
    lat_end
        起始经度，可参考 arcmap 生成的栅格文件 source 属性
    lon_start
         纬度同上
    lon_end
         纬度同上
    degree
        输出栅格的一个网格的度数
    output_file
        输出 .tif 文件的路径

    Returns
    -------
    None
    """

    nx, ny = array.shape
    mag_grid = np.reshape(array, (nx, ny), order="F")  # !!!
    mag_grid = np.float64(mag_grid)
    lats = np.linspace(start=lat_start, stop=lat_end, num=mag_grid.shape[0])
    lons = np.linspace(start=lon_start, stop=lon_end, num=mag_grid.shape[1])
    assert len(lats) == mag_grid.shape[0]
    assert len(lons) == mag_grid.shape[1]
    xres = lons[1] - lons[0]
    yres = lats[1] - lats[0]
    ysize = len(lats)
    xsize = len(lons)
    driver = gdal.GetDriverByName("GTiff")
    ds = driver.Create(output_file, xsize, ysize, 1, gdal.GDT_Float32)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    gt = [lon_start, xres, 0, lat_start, 0, yres]
    ds.SetGeoTransform(gt)
    outband = ds.GetRasterBand(1)
    outband.SetStatistics(
        np.min(mag_grid), np.max(mag_grid), np.average(mag_grid), np.std(mag_grid)
    )
    outband.WriteArray(mag_grid)
    ds = None


def shp_id(shpfile: str):
    """
    Get basin id from shapefile name

    Parameters
    ----------
    shpfile
        shapefile path e.g. ./0000.shp
    Returns
    -------
    str
        shapefile id e.g. 0000
    """
    # original code
    #return re.findall(r"[\d]+", shpfile)[-1]

    # adaptive code
    # Remove the file extension ".shp" and then truncate what comes after "basin_"
    basename = os.path.splitext(shpfile)[0]
    result = basename.split('basin_', 1)[-1]
    return result


def absolute_file_paths(directory):
    """
    List all files in the directory and its all subdirectories

    Parameters
    ----------
    directory
        the dir path

    Returns
    -------
    list
        all files in the directory and its all subdirectories
    """

    def nest(nest_directory):
        for path, _, filenames in os.walk(nest_directory):
            for f in filenames:
                yield os.path.abspath(os.path.join(path, f))

    return list(nest(directory))


def reproject_tif(src_tif: str, out_tif: str, out_crc="EPSG:4326"):
    """
    Reproject src_tif file to out_tif with out_crc

    Parameters
    ----------
    src_tif
        source tif file path
    out_tif
        output tif file path
    out_crc
        target coordination system; default is EPSG:4326

    Returns
    ----------
    None
    """
    with rasterio.open(src_tif) as src:
        transform, width, height = calculate_default_transform(
            src.crs, out_crc, src.width, src.height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update(
            {"crs": out_crc, "transform": transform, "width": width, "height": height}
        )

        with rasterio.open(out_tif, "w", **kwargs) as dst:
            for i in range(1, src.count + 1):
                reproject(
                    source=rasterio.band(src, i),
                    destination=rasterio.band(dst, i),
                    src_transform=src.transform,
                    src_crs=src.crs,
                    dst_transform=transform,
                    dst_crs=out_crc,
                    resampling=Resampling.nearest,
                )


def merge_tifs(tif_files: list, outfile: str):
    """
    Merge multiple tif files

    Parameters
    ----------
    tif_files
        list of .tif file paths
    outfile
        output file path

    Returns
    -------
    None
    """
    src_files_to_mosaic = []
    for fp in tif_files:
        src = rasterio.open(fp)
        src_files_to_mosaic.append(src)
    mosaic, out_trans = merge(src_files_to_mosaic)
    out_meta = src.meta.copy()
    out_meta.update(
        {
            "driver": "GTiff",
            "height": mosaic.shape[1],
            "width": mosaic.shape[2],
            "transform": out_trans,
            "crs": "EPSG:4326",
        }
    )
    with rasterio.open(outfile, "w", **out_meta) as dest:
        dest.write(mosaic)


def extract_raster_by_shape_file(
    raster: str, shape_file: str, output_file=None, nodata=-9999
) -> np.array:
    """
    Extract data in the given shapefile from a raster file

    Parameters
    ----------
    raster
        EPSG:4326 .tif file path
    shape_file
        EPSG:4326 .shp file path
    output_file
        output .tif file path; default is None which means no output
    nodata
        nodata's value; default is -9999

    Returns
    -------
    np.array
        the raster values array on the given polygon of the shapefile
    """
    with fiona.open(shape_file, "r") as shapefile:
        shapes = [feature["geometry"] for feature in shapefile]
    with rasterio.open(raster) as src:
        out_image, out_transform = rasterio.mask.mask(
            src, shapes, nodata=nodata, crop=True
        )
        out_meta = src.meta
    if output_file is None:
        return out_image
    else:
        out_meta.update(
            {
                "driver": "GTiff",
                "height": out_image.shape[1],
                "width": out_image.shape[2],
                "transform": out_transform,
            }
        )
        with rasterio.open(output_file, "w", **out_meta) as dest:
            dest.write(out_image)
        return out_image


def zonal_stats_singletif(
    tif_file: str, shape_file: str, valid_min=None, valid_max=None
) -> float:
    """
    Input one .shp file and one .tif file, and make zonal statistics

    Parameters
    ----------
    tif_file
        the input tif file
    shape_file
        the input shapefile

    Returns
    -------
    float
        zonal statistics
    """
    res = extract_raster_by_shape_file(tif_file, shape_file).flatten()
    res = res[res != -9999]
    res = res[~np.isnan(res)]
    if valid_min is not None:
        res = res[res > valid_min]
    if valid_max is not None:
        res = res[res < valid_max]
    if len(res) > 0:
        return np.mean(res)
    else:
        return np.nan
