import datetime
import os.path
import shutil
import subprocess
import sys

import pandas as pd
from tqdm import tqdm

from utils import *

"""
According to MODIS product，calculate NDVI/LAI basin mean value daily series

reference:
https://lpdaac.usgs.gov/products/mcd15a3hv006/
https://lpdaac.usgs.gov/products/mod13q1v006/

Requirement:
(1) MODIS data
├── MOD13Q1/MCD15A3H
|   ├── MCD15A3H.A2002185.h22v04.006.2015149102803.hdf
|   ├── MCD15A3H.A2002186.h22v04.006.2015149102803.hdf
|   ├── MCD15A3H.A2002187.h22v04.006.2015149102803.hdf
|   ├── ......
(2) basins shapefiles
├── folder_shp
|   ├── basin_xxxx.shp
|   ├── ...
"""


def get_qualified_hdf_files_from_folder(folder: str, modis_product: str, zones: list):
    """

    :param folder: modis hdfs folder
    :param modis_product: modis product id: 'MCD15A3H' for LAI; 'MOD13Q1' for NDVI
    :param zones: modis tiles
                  https://lpdaac.usgs.gov/data/get-started-data/collection-overview/missions/modis-overview/
                  e.g. [h25v06, h25v07]
    :return: qualified hdfs list
    """
    if zones == "all":
        res = [file for file in absolute_file_paths(folder) if file.endswith(".hdf")]
    else:
        res = [
            file
            for file in absolute_file_paths(folder)
            if file.endswith(".hdf")
            and (get_info_from_modis_hdf(file)["zones"] in zones)
        ]
    for file in res:
        assert os.path.basename(file).split(".")[0] == modis_product
    return res


def get_info_from_modis_hdf(file_path) -> dict:
    """
    file_path: 'MCD15A3H.A2018017.h25v06.006.2018023210623.hdf' or absolute path

    > get_info_from_modis_hdf('MCD15A3H.A2018017.h25v06.006.2018023210623.hdf')
    """
    res = {}
    date = os.path.basename(file_path).split(".")[1]
    year = int(date[1:-3])
    day_of_year = int(date[-3:])
    res["date"] = datetime.datetime(year, 1, 1) + datetime.timedelta(day_of_year - 1)
    res["product"] = os.path.basename(file_path).split(".")[0]
    res["zones"] = os.path.basename(file_path).split(".")[2]
    return res


def get_info_from_modis_tif(file_path) -> dict:
    """
    Get inforamtion of MODIS tif file

    Parameters
    ----------
    file_path
        "MCD12Q1.A2018001.h25v04.006.2019200013451_08.tif" or absolute path

    Returns
    -------
    dict
        date, product, zones, feature
    """
    print(file_path)
    res = {}
    date = os.path.basename(file_path).split(".")[1]
    print(date)
    year = int(date[1:-3])
    day_of_year = int(date[-3:])
    res["date"] = datetime.datetime(year, 1, 1) + datetime.timedelta(day_of_year - 1)
    res["product"] = os.path.basename(file_path).split(".")[0]
    res["zones"] = os.path.basename(file_path).split(".")[2]
    res["feature"] = re.findall("_(\d+)", file_path)[-1]
    return res


def group_hdf_files_by_date(hdf_files: list):
    """

    :param hdf_files: hdf files
    :return: group hdf files by dare
    """
    dates = [get_info_from_modis_hdf(file)["date"] for file in hdf_files]
    unique_dates = np.unique(dates)
    res = {}
    for date in unique_dates:
        res[date] = [
            file for file in hdf_files if date == get_info_from_modis_hdf(file)["date"]
        ]
    return res


def hdf_to_tif(hdf_file: str, output_dir: str):
    """
    Transform hdf file to tif file

    Parameters
    ----------
    hdf_file
        hdf file
    output_dir
        convert hdf file to tifs
    Returns
    -------
    str
        tif file's name
    """
    cwd = os.getcwd()
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    os.chdir(output_dir)
    name = os.path.basename(hdf_file)[:-4] + ".tif"
    command = f"gdal_translate -sds -of GTiff {hdf_file} {name}"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.chdir(cwd)
    return name


def gdal_downsample_tif(tif_file: str, working_dir: str, percent: int):
    """
    Downsample tif

    Parameters
    ----------
    tif_file
        tif file path
    working_dir
        data processing root dir
    percent
        downsample percent
    Returns
    -------
    None
    """
    cwd = os.getcwd()
    os.chdir(working_dir)
    ori_name = os.path.basename(tif_file)
    new_name = os.path.basename(ori_name)[:-4] + "_tmp.tif"
    os.rename(tif_file, new_name)
    command = f"gdal_translate -outsize {percent}% GTiff {new_name} {ori_name}"
    process = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    os.remove(new_name)
    os.chdir(cwd)


def group_tif_files_by_date_feature(files: list):
    """
    Group tifs by date and feature (there are some features in the original MODIS HDF file)

    Parameters
    ----------
    files
        tif files

    Returns
    -------
    dict
        group tifs by date and feature
    """
    dates = [get_info_from_modis_tif(file)["date"] for file in files]
    unique_dates = np.unique(dates)
    features = [get_info_from_modis_tif(file)["feature"] for file in files]
    unique_features = np.unique(features)
    res = {}
    for date in unique_dates:
        res[date] = {}
        for feature in unique_features:
            res[date][feature] = [
                file
                for file in files
                if date == get_info_from_modis_tif(file)["date"]
                and feature == get_info_from_modis_tif(file)["feature"]
            ]
    return res, unique_features


def get_84_tifs(folder: str) -> list:
    """
    Return only tifs in wgs84

    Parameters
    ----------
    folder
        tifs folder

    Returns
    -------
    list
        tifs in wgs84
    """
    files = absolute_file_paths(folder)
    return [file for file in files if "tif84" in file]


def zonal_stats(tif_file: str, shape_file: str, valid_min, valid_max) -> dict:
    """
    Zonal statistics

    Parameters
    ----------
    tif_file
        tif file path
    shape_file
        shp file path
    valid_min
        NDVI: [-2000, 10000]; LAI: [0, 100]
    valid_max
        NDVI: [-2000, 10000]; LAI: [0, 100]

    Returns
    -------
    dict
        {'mean': xx, 'max': xx, 'min': xx}
    """
    res = extract_raster_by_shape_file(tif_file, shape_file).flatten()
    res[res > valid_max] = valid_max + 1
    res[res < valid_min] = valid_max + 1
    res = res[res != valid_max + 1]
    res = res[~np.isnan(res)]
    if len(res) > 0:
        return {"mean": np.mean(res), "max": np.max(res), "min": np.min(res)}
    else:
        return {"mean": 0, "max": 0, "min": 0}


def clear_dir(folder: str):
    """
    Clear all cache directory

    Parameters
    ----------
    folder
        folder path
    Returns
    -------
    None
    """
    shutil.rmtree(folder)
    print(f"{folder} cleared")


class Modis:
    def __init__(self, hdf_folder, tmp_folder, product, zones):
        """
        Process MODIS product data

        Parameters
        ----------
        hdf_folder
            directory where we put hdf files
        tmp_folder
            directory where we put tif files
            after downsampling the tif files, they will be removed, hence, we call the dir "tmp_folder"
        product
            MODIS product name
        zones
            MODIS tiles' zones
        """
        self.hdf_folder = hdf_folder
        self.tmp_folder = tmp_folder
        self.product = product
        self.zones = zones
        self.merged_tif_names = []
        print("hdf_folder: ", hdf_folder)
        print("tmp_folder: ", tmp_folder)
        print("product: ", product)
        print("zones: ", zones)
        print("---------------------")

    def get_merged_tifs(
        self,
        merged_tifs_folder: str,
        feature_name: str,
        feature_index: str,
        downsampled_cache: str,
    ) -> pd.DataFrame:
        """
        Merge tif files to one

        Parameters
        ----------
        merged_tifs_folder
            the merged tif files' directory;
            after downsampling, the files in tif_cache will be removed,
            hence we directly use tif_cache as merged_tifs_folder
        feature_name
            LAI
        feature_index
            there are 6 features in the MCD15A3H.006 HDF file, the second one is LAI (1 in HDF / 2 in our directory)
        downsampled_cache
            downsampled tif files' directory

        Returns
        -------
        pd.DataFrame
            merged tif data
        """
        print("merged_tifs_folder: ", merged_tifs_folder)
        print("feature_index: ", feature_index)
        print("feature_name: ", feature_name)
        print("----------------------")
        print("clear the tmp folder, if not exists, creat one")
        files = absolute_file_paths(self.hdf_folder)
        # files = get_qualified_hdf_files_from_folder(self.hdf_folder, self.product, self.zones)
        print("convert hdf to tifs")
        for file in tqdm(files, position=0, leave=True, file=sys.stdout):
            hdf_to_tif(file, self.tmp_folder)
        print("resize tifs to 50%")
        for file in tqdm(absolute_file_paths(self.tmp_folder)):
            if file.endswith(".tif"):
                gdal_downsample_tif(file, downsampled_cache, 50)
        print("reproject tifs to WGS:84")
        tifs = [
            file
            for file in absolute_file_paths(downsampled_cache)
            if file.endswith(".tif")
        ]
        for file in tqdm(tifs, position=0, leave=True, file=sys.stdout):
            name = file + "84.tif"
            reproject_tif(file, name)
        tifs_84 = get_84_tifs(downsampled_cache)
        groups, unique_features = group_tif_files_by_date_feature(tifs_84)
        print("merge tifs and stats")
        if not os.path.isdir(merged_tifs_folder):
            os.makedirs(merged_tifs_folder)
        merged_tifs = {}
        for date in tqdm(groups.keys(), position=0, leave=True, file=sys.stdout):
            merged_tifs[date] = {}
            merged_tif_name = f"{self.product}-{date.year}.{date.month}.{date.day}-{feature_name}-merged.tif"
            merged_tif_name = os.path.join(merged_tifs_folder, merged_tif_name)
            merge_tifs(groups[date][feature_index], merged_tif_name)
            merged_tifs[date][feature_name] = merged_tif_name
            self.merged_tif_names.append(merged_tif_name)
        return pd.DataFrame(merged_tifs).T

    def clear_tmp(self):
        shutil.rmtree(self.tmp_folder)


def get_hdf_product(file):
    """

    :param file: e.g. ./MCD15A3H.A2002185.h23v03.006.2015149105852.hdf
    :return: 'MCD15A3H'
    """
    return os.path.basename(file).split(".")[0]


def get_hdf_date(file):
    """
    get the date of the hdf file

    Parameters
    ----------
    file
        e.g. ./MCD15A3H.A2002185.h23v03.006.2015149105852.hdf

    Returns
    -------
    datetime.datetime
        the date of the hdf file; datetime.datetime(2002, M, D)
    """
    tmp = os.path.basename(file).split(".")[1]
    year = int(tmp[1:5])
    days = int(tmp[5:])
    return datetime.datetime(year, 1, 1) + datetime.timedelta(days - 1)


def summary_year(
    year, data_root, out_dir, hdf_dir, shp_dir, tif_dir, downsampled_cache, modis_dict
):
    """
    Summary for one year

    Parameters
    ----------
    year
        specify the year to calculate
    data_root
        modis lai/ndvi data root dir, e.g. ./MOD13Q1
    out_dir
        output dir, e.g. ./data
    hdf_dir
        temporary hdf files' directory
    shp_dir
        directory where we put basins' shpefiles
    tif_dir
        directory where we put intermediate tif files
    downsampled_cache
        downscaled hdf files' directory
    modis_dict
        parameters dict of LAI or NDVI:
        LAI -- {"product": "MCD15A3H", "feature_name": "LAI", "feature_index": "2", "valid_min": 0, "valid_max": 100}
        NDVI -- {"product": "MOD13Q1", "feature_name": "NDVI", "feature_index": "1", "valid_min": -2000, "valid_max": 10000}

    Returns
    -------
    None
    """
    start = datetime.datetime(year, 1, 1)
    end = datetime.datetime(year, 12, 31)

    if os.path.exists(hdf_dir):
        clear_dir(hdf_dir)
        os.makedirs(hdf_dir)
    else:
        os.makedirs(hdf_dir)

    if os.path.exists(tif_dir):
        clear_dir(tif_dir)
        os.makedirs(tif_dir)
    else:
        os.makedirs(tif_dir)

    if os.path.exists(downsampled_cache):
        clear_dir(downsampled_cache)
        os.makedirs(downsampled_cache)
    else:
        os.makedirs(downsampled_cache)

    files = absolute_file_paths(data_root)
    files = [file for file in files if file.endswith(".hdf")]

    res = {}
    hdf_files_dates = np.unique([get_hdf_date(file) for file in files])
    for date in tqdm(hdf_files_dates):
        if not start <= date <= end:
            continue

        print(date, "...")
        for file in files:
            if get_hdf_date(file) == date:
                shutil.copyfile(file, os.path.join(hdf_dir, os.path.basename(file)))

        modis = Modis(
            hdf_folder=hdf_dir,
            tmp_folder=tif_dir,
            product=modis_dict["product"],
            zones="all",
        )
        modis.get_merged_tifs(
            merged_tifs_folder=tif_dir,
            feature_name=modis_dict["feature_name"],
            feature_index=modis_dict["feature_index"],
            downsampled_cache=downsampled_cache,
        )

        for shape_file in [
            file for file in absolute_file_paths(shp_dir) if file.endswith(".shp")
        ]:
            shape_id = shp_id(shape_file)
            tmp_res = zonal_stats(
                tif_file=modis.merged_tif_names[0],
                shape_file=shape_file,
                valid_min=modis_dict["valid_min"],
                valid_max=modis_dict["valid_max"],
            )
            if shape_id not in res:
                res[shape_id] = {}
            res[shape_id][date] = tmp_res["mean"]

        clear_dir(hdf_dir)
        os.makedirs(hdf_dir)
        clear_dir(tif_dir)
        os.makedirs(tif_dir)
        clear_dir(downsampled_cache)
        os.makedirs(downsampled_cache)

    for key in res:
        year_dir = os.path.join(out_dir, str(year))
        if not os.path.isdir(year_dir):
            os.makedirs(year_dir)
        pd.DataFrame(res[key], index=[0]).T.to_csv(
            os.path.join(year_dir, f"{key}.csv"), header=None
        )
