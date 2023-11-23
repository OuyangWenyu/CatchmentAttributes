"""Here are CAMELS-CC's scripts"""
import argparse
import os
import sys
import pandas as pd
import geopandas as gpd
from tqdm import tqdm

sys.path.append("..")

import definitions
from basin_era5_process import trans_era5_land_to_camels_format
from climate import (
    series_mean,
    high_prec_freq,
    high_prec_dur,
    high_prec_timing,
    low_prec_dur,
    low_prec_freq,
    low_prec_timing,
    frac_snow_daily, p_seasonality,
)
from glim import Glim
from igbp import igbp_stats
from rooting_depth import DepthMapper, root_depth_50_99_stats
from topo_elev import elev_mean, slope_mean, merge_and_reproject_dems
from topo_shape import basin_topo_stats
from modis import summary_year
from utils import absolute_file_paths, reproject_tif, shp_id, zonal_stats_singletif
from soil import binary2tif, tif_from_nc, all_soil_depth_mean_weight_in_soilgrids250
from permeability_porosity import GLHYMPS


def res_to_df(res):
    """
    Turn the result dict to a DataFrame with sorted gage_ids

    Parameters
    ----------
    res
        result attributes

    Returns
    -------
    pd.DataFrame
        the attributes dataFrame to be saved
    """
    df_res = pd.DataFrame(res).T
    df_res["gage_id"] = df_res.index
    df_res_order = df_res[["gage_id"] + df_res.columns[:-1].values.tolist()]
    df_res_order_sort = df_res_order.sort_values(by=["gage_id"])
    return df_res_order_sort


def basin_mean_forcing(args):
    """
    Main function which is called from the command line. Entrypoint for training all ML models.
    """
    print(
        "------------------------Calculate basin mean forcing from source data---------------------------------"
    )
    input_dir = os.path.join(definitions.DATASET_DIR, "ERA5_LAND")
    if not os.path.isdir(input_dir):
        raise NotADirectoryError("Please put your forcing data to data/ERA5_LAND")
    output_dir = os.path.join(definitions.DATASET_DIR, "basin_mean_forcing")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    region = args.region
    years = list(range(int(args.year_range[0]), int(args.year_range[1])))
    gage_dict_file = os.path.join(definitions.DATASET_DIR, region + "_name.txt")
    if not os.path.isfile(gage_dict_file):
        raise FileNotFoundError(
            "No such id file, please set " + gage_dict_file + ". A template could be seen in data/camels_mr_name.txt")
    gage_dict = pd.read_csv(gage_dict_file, dtype={"gage_id": str}).to_dict(orient="list")
    for year in years:
        trans_era5_land_to_camels_format(input_dir, output_dir, gage_dict, region, year)
    print("------------------------Finished---------------------------------")


def climate_app():
    print(
        "------------------------Calculate basin climate attributes---------------------------------"
    )
    forcing_dir = os.path.join(definitions.DATASET_DIR, "basin_mean_forcing")
    if not os.path.isdir(forcing_dir):
        raise NotADirectoryError("CAMELS-CC need your forcing data!")
    output_dir = os.path.join(definitions.DATASET_DIR, "attribute")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    files = [x for x in absolute_file_paths(forcing_dir)]
    res = {}
    #Do this for every txt file
    for file in tqdm(files):
        #This place intercepts a _ but the new data gage_id has a - inside. this is very troublesome need to change this place
        filename = os.path.basename(file)
        target_extension = "_lump_era5_land_forcing"
        # Find the target extension at the beginning of the filename
        start_index = filename.rfind(target_extension)
        # Use a slice to capture the part of the filename after the last underscore
        name = filename[:start_index]

        #The following line of code is original
        #name = file.split(os.path.sep)[-1].split("_")[0]

        df = pd.read_csv(file, sep="\s+")
        pre = df["total_precipitation"]
        # evp data in ERA5-LAND are negative
        pet = df["potential_evaporation"] * (-1)
        # there are some abnormal values
        aet = df[(df["total_evaporation"] > -1) & (df["total_evaporation"] < 0)]["total_evaporation"] * (-1)
        # m/day -> mm/day
        pet_mean = series_mean(pet) * 1000
        p_mean = series_mean(pre) * 1000
        p_s = p_seasonality(df)[0]
        res[name] = {
            "p_mean": p_mean,
            "pet_mean": pet_mean,
            "aet_mean": series_mean(aet) * 1000,
            "aridity": pet_mean / p_mean,
            "p_seasonality": p_s,
            "high_prec_freq": high_prec_freq(pre),
            "high_prec_dur": high_prec_dur(pre),
            "high_prec_timing": high_prec_timing(df),
            "low_prec_freq": low_prec_freq(pre),
            "low_prec_dur": low_prec_dur(pre),
            "low_prec_timing": low_prec_timing(df),
            "frac_snow_daily": frac_snow_daily(df),
        }
        # p_seasonality's calculation is too slow, so we quit it now
    df_res_order_sort = res_to_df(res)
    df_res_order_sort.to_csv(os.path.join(output_dir, "climate.csv"), index=None,mode='a')
    print("------------------------Finished---------------------------------")


def glim_app():
    print(
        "------------------------Calculate basin geology attributes---------------------------------"
    )
    path_glim_tif = os.path.join(
        definitions.DATASET_DIR, "GLIM", "glim_clip_raster.tif"
    )
    glim_raster_tif = os.path.join(definitions.DATASET_DIR, "GLIM", "GlimRaster.tif")
    if not os.path.isfile(glim_raster_tif):
        reproject_tif(path_glim_tif, glim_raster_tif, out_crc="EPSG:4326")

    glim_cate_number_mapping_file = os.path.join(
        definitions.DATASET_DIR, "GLIM", "GLiMCateNumberMapping.csv"
    )
    if not os.path.isfile(glim_cate_number_mapping_file):
        litho_index_field_file = os.path.join(
            definitions.DATASET_DIR, "GLIM", "xx_index_field.shp"
        )
        if not os.path.isfile(litho_index_field_file):
            raise FileNotFoundError(
                "Please see README.md to make sure you have generated xx_index_field.shp"
            )
        gdf = gpd.read_file(litho_index_field_file)
        gdf[["OBJECTID", "xxValue", "xx"]].to_csv(
            glim_cate_number_mapping_file, index=None
        )

    short2long_name_txt = os.path.join(
        definitions.DATASET_DIR, "glim_name_short_long.txt"
    )

    glimer = Glim(
        glim_raster_tif=glim_raster_tif,
        glim_cate_number_mapping_file=glim_cate_number_mapping_file,
        short2long_name_txt=short2long_name_txt,
    )

    res = {}
    basin_shape_files_dir = os.path.join(definitions.DATASET_DIR, "shapefiles")
    for shape_file in tqdm(
            [
                file
                for file in absolute_file_paths(basin_shape_files_dir)
                if file.endswith(".shp")
            ]
    ):
        res[shp_id(shape_file)] = glimer.extract_basin_attributes_glim(
            shape_file=shape_file
        )
    output_dir = os.path.join(definitions.DATASET_DIR, "attribute")
    df_res_order_sort = res_to_df(res)
    df_res_order_sort.to_csv(os.path.join(output_dir, "geology.csv"), index=None)
    print("------------------------Finished---------------------------------")


def permeability_porosity_app():
    print(
        "----------------------Calculate basin permeability and porosity-------------------------------"
    )
    permeability_no_permafrost_raster_tif = os.path.join(
        definitions.DATASET_DIR, "9_code_data", "processed_permeability.tif"
    )
    porosity_raster_tif = os.path.join(
        definitions.DATASET_DIR, "9_code_data", "processed_porosity.tif"
    )
    nan_value = 65535

    glhympser = GLHYMPS(
        permeability_no_permafrost_raster_tif, porosity_raster_tif, nan_value=nan_value
    )

    res = {}
    basin_shape_files_dir = os.path.join(definitions.DATASET_DIR, "shapefiles")
    for shape_file in tqdm(
            [
                file
                for file in absolute_file_paths(basin_shape_files_dir)
                if file.endswith(".shp")
            ]
    ):
        res[shp_id(shape_file)] = glhympser.zonal_stats_glhymps(shape_file=shape_file)
    res = pd.DataFrame(res).T.reset_index().rename(columns={"index": "gage_id"})

    output_dir = os.path.join(definitions.DATASET_DIR, "attribute")
    df_res_order_sort = pd.DataFrame(res).sort_values(by=["gage_id"])
    df_res_order_sort.to_csv(
        os.path.join(output_dir, "permeability_porosity.csv"), index=None
    )
    print("------------------------Finished---------------------------------")


def igbp_app():
    print(
        "------------------------Calculate basin land cover attributes---------------------------------"
    )
    igbp_tif = os.path.join(definitions.DATASET_DIR, "processed_igbp.tif")
    shp_dir = os.path.join(definitions.DATASET_DIR, "shapefiles")
    out = os.path.join(definitions.DATASET_DIR, "attribute", "land_cover.csv")

    res = {}
    for shape_file in tqdm(
            file for file in absolute_file_paths(shp_dir) if file.endswith(".shp")
    ):
        res[shp_id(shape_file)] = igbp_stats(shapefile=shape_file, igbp_tif=igbp_tif)
    df_res_order_sort = res_to_df(res)
    df_res_order_sort.to_csv(out, index=None)
    print("------------------------Finished---------------------------------")


def root_depth_app():
    print(
        "------------------------Calculate basin root depth attributes---------------------------------"
    )
    igbp_tif = os.path.join(definitions.DATASET_DIR, "processed_igbp.tif")
    root_depth = os.path.join(definitions.DATASET_DIR, "calculated_root_depth.txt")
    shp_dir = os.path.join(definitions.DATASET_DIR, "shapefiles")
    out = os.path.join(definitions.DATASET_DIR, "attribute", "root_depth.csv")
    depth_mapper = DepthMapper(root_depth)

    res = {}
    for shape_file in tqdm(
            file for file in absolute_file_paths(shp_dir) if file.endswith(".shp")
    ):
        res[shp_id(shape_file)] = root_depth_50_99_stats(
            shape_file, igbp_tif, depth_mapper
        )
    df_res_order_sort = res_to_df(res)
    df_res_order_sort.to_csv(out, index=None)
    print("------------------------Finished---------------------------------")


def topo_elev_app():
    print(
        "------------------------Calculate basin topography elevation and slope---------------------------------"
    )
    dem_folder = os.path.join(definitions.DATASET_DIR, "DEM")
    shp_folfer = os.path.join(definitions.DATASET_DIR, "shapefiles")
    outpath = os.path.join(definitions.DATASET_DIR, "attribute", "topo_elev_slope.csv")
    topo_dir = os.path.join(definitions.DATASET_DIR, "TOPO")
    if not os.path.isdir(topo_dir):
        os.makedirs(topo_dir)
    tmp_merged = os.path.join(topo_dir, "merge_cache.tif")
    tmp_reprojected = os.path.join(topo_dir, "reproject_cache.tif")
    tmp_slope = os.path.join(topo_dir, "slope_cache.tif")
    res = []
    shps = [file for file in absolute_file_paths(shp_folfer) if file.endswith(".shp")]
    for shpfile in tqdm(shps):
        merge_and_reproject_dems(shpfile, dem_folder, tmp_merged, tmp_reprojected)
        tmp_res = {
            "gage_id": shp_id(shpfile),
            "elev": elev_mean(shpfile, tmp_reprojected),
            "slope": slope_mean(shpfile, tmp_reprojected, tmp_slope),
        }
        res.append(tmp_res)
    df_res_order_sort = pd.DataFrame(res).sort_values(by=["gage_id"])
    df_res_order_sort.to_csv(outpath, index=None)
    print("------------------------Finished---------------------------------")


def topo_shape_app():
    print(
        "------------------------Calculate basin topography shape attributes---------------------------------"
    )
    shp_folfer = os.path.join(definitions.DATASET_DIR, "shapefiles")
    stream_shps = os.path.join(definitions.DATASET_DIR, "TOPO", "as_streams.shp")
    out_path = os.path.join(
        definitions.DATASET_DIR, "attribute", "topo_shape_factors.csv"
    )

    basin_shps = [file for file in absolute_file_paths(shp_folfer) if ".shp" in file]
    topo_stats = basin_topo_stats(basin_shps=basin_shps, stream_shps=stream_shps)
    df_res_order_sort = res_to_df(topo_stats)
    df_res_order_sort.to_csv(out_path, index=None)
    print("------------------------Finished---------------------------------")


def modis_app(args):
    print(
        "------------------------Calculate basin mean LAI or max NDVI---------------------------------"
    )
    modis_type = args.modis_type
    year_range = [int(args.year_range[0]), int(args.year_range[1])]
    if modis_type == "lai":
        modis_dict = {
            "product": "MCD15A3H",
            "feature_name": "LAI",
            "feature_index": "2",
            "valid_min": 0,
            "valid_max": 100,
        }
    elif modis_type == "ndvi":
        modis_dict = {
            "product": "MOD13Q1",
            "feature_name": "NDVI",
            "feature_index": "01",
            "valid_min": -2000,
            "valid_max": 10000,
        }
    else:
        raise NotImplementedError("We didn't provide tools for this MODIS product!!!")
    hdf_dir = os.path.join(
        definitions.DATASET_DIR, modis_dict["feature_name"], "hdf_cache"
    )
    tif_dir = os.path.join(
        definitions.DATASET_DIR, modis_dict["feature_name"], "tif_cache"
    )
    downsampled_cache = os.path.join(
        definitions.DATASET_DIR, modis_dict["feature_name"], "downsampled_cache"
    )
    shp_dir = os.path.join(definitions.DATASET_DIR, "shapefiles")
    data_root = os.path.join(definitions.DATASET_DIR, modis_dict["product"])
    out_dir = os.path.join(definitions.DATASET_DIR, modis_dict["feature_name"])
    for year in range(year_range[0], year_range[1]):
        summary_year(
            year,
            data_root,
            out_dir,
            hdf_dir,
            shp_dir,
            tif_dir,
            downsampled_cache,
            modis_dict=modis_dict,
        )
    print("------------------------Finished---------------------------------")


def soil_app():
    print(
        "----------------------Calculate basin soil attribute-------------------------------"
    )

    # Here we only use tif files from https://files.isric.org/soilgrids/former/2017-03-10/data/
    # we mainly process these attributes:
    # bdticm: Depth to bedrock
    # sand: Soil sand content (%), average over all layers
    # clay: Soil clay content (%), average over all layers
    soil_source_dir = os.path.join(definitions.DATASET_DIR, "soil_source_data")
    # zonal stats
    res = {}
    files = [x for x in absolute_file_paths(soil_source_dir) if x.endswith(".tif")]
    shapefiles_dir = os.path.join(definitions.DATASET_DIR, "shapefiles")
    shps = [x for x in absolute_file_paths(shapefiles_dir) if x.endswith(".shp")]
    for file in files:
        if not ("_downscaled" in file):
            continue
        try:
            for shp in shps:
                if not shp_id(shp) in res:
                    res[shp_id(shp)] = {}
                var_name = (
                    os.path.basename(file).split(".")[0].replace("_downscaled", "")
                )
                res[shp_id(shp)][var_name] = zonal_stats_singletif(
                    file, shp, valid_min=0, valid_max=None
                )
        except Exception as e:
            print(e)
            continue
    res = pd.DataFrame(res).T.reset_index().rename(columns={"index": "gage_id"})
    # most attributes in SoilGrids250m include values of 7 standard depths, we calculate mean value
    res_new = all_soil_depth_mean_weight_in_soilgrids250(res)
    output_dir = os.path.join(definitions.DATASET_DIR, "attribute")
    df_res_order_sort = pd.DataFrame(res_new).sort_values(by=["gage_id"])
    df_res_order_sort.to_csv(os.path.join(output_dir, "soil.csv"), index=None)
    print("------------------------Finished---------------------------------")


# python app.py --catch_attr basin_mean_forcing --year_range 2014 2022 --region camels_mr
# python app.py --catch_attr climate
# python app.py --catch_attr geology
# python app.py --catch_attr permeability_porosity
# python app.py --catch_attr land_cover
# python app.py --catch_attr root_depth
# python app.py --catch_attr topo_elev
# python app.py --catch_attr topo_shape
# python app.py --catch_attr modis --modis_type lai --year_range 2015 2016
# python app.py --catch_attr modis --modis_type ndvi --year_range 2015 2016
# python app.py --catch_attr soil
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Calculate Attributes of Catchments"
    )
    parser.add_argument(
        "--catch_attr",
        dest="catch_attr",
        help="Chose one function to run",
        default="soil",
        type=str,
    )
    parser.add_argument(
        "--modis_type",
        dest="modis_type",
        help="Chose one type of MODIS, here we provide lai and ndvi",
        default="ndvi",
        type=str,
    )
    parser.add_argument(
        "--year_range",
        dest="year_range",
        help="The start and end years (right open interval)",
        default=[2015, 2016],
        nargs="+",
    )
    parser.add_argument(
        "--region",
        dest="region",
        help="The name of region including all target basins",
        default="camels_mr",
        type=str,
    )
    the_args = parser.parse_args()
    if the_args.catch_attr == "basin_mean_forcing":
        basin_mean_forcing(the_args)
    elif the_args.catch_attr == "climate":
        climate_app()
    elif the_args.catch_attr == "geology":
        glim_app()
    elif the_args.catch_attr == "permeability_porosity":
        permeability_porosity_app()
    elif the_args.catch_attr == "land_cover":
        igbp_app()
    elif the_args.catch_attr == "root_depth":
        root_depth_app()
    elif the_args.catch_attr == "topo_elev":
        topo_elev_app()
    elif the_args.catch_attr == "topo_shape":
        topo_shape_app()
    elif the_args.catch_attr == "modis":
        modis_app(the_args)
    elif the_args.catch_attr == "soil":
        soil_app()
    else:
        raise NotImplementedError("We have not provided this function yet!!!")