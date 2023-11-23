"""
Process ERA5 land data for basins
"""
import datetime
import fnmatch
import os
from typing import Union

import numpy as np
import pandas as pd
from dateutil import tz


def utc_to_local(
    utc_time: Union[str, np.datetime64],
    local_tz: str = "Asia/Hong_Kong",
    time_format="%Y-%m-%dT%H:%M:%S",
):
    """
    Transform UTC time to local time.

    Code comes from: https://stackoverflow.com/questions/4770297/convert-utc-datetime-string-to-local-datetime

    Parameters
    ----------
    utc_time
        the UTC time
    local_tz
        the local time zone; default is "Asia/Hong_Kong"
    time_format
        the format of utc_time_str and return value; default is "%Y-%m-%dT%H:%M:%S"

    Returns
    -------
    str
        return local time string; its format is time_format
    """
    from_zone = tz.gettz("UTC")
    to_zone = tz.gettz(local_tz)
    if type(utc_time) is np.datetime64:
        utc_time = np.datetime_as_string(utc_time)
    utc = datetime.datetime.strptime(utc_time, time_format)
    # Tell the datetime object that it's in UTC time zone since
    # datetime objects are 'naive' by default
    utc = utc.replace(tzinfo=from_zone)
    # Convert time zone
    central = utc.astimezone(to_zone)
    result_str = central.strftime(time_format)
    return result_str


def trans_era5_land_to_camels_format(
    era5_land_dir, output_dir, gage_dict, region, year, time_zone="Asia/Hong_Kong"
):
    """
    Transform hourly forcing data of ERA5-LAND downloaded from GEE to the format in CAMELS.

    ERA5-LAND's forcing data is hourly, so we resample it to daily in GEE,
    In addition, we will use them for basins in China, so we chose the time zone "Asia/Hong_Kong"
    The GEE code: https://code.earthengine.google.com/b40bfe7529ec9df928f5abfc029f4e4e

    Parameters
    ----------
    era5_land_dir
        the original data's directory
    output_dir
        the transformed data's directory
    gage_dict
        a dict containing gage's ids and the correspond HUC02 ids
    region
        we named the file downloaded from GEE as era5_land_<region>_mean_<year>.csv,
        because we use GEE code to generate data for each year for each shape file (region) containing some basins.
        For example, if we use the basins' shpfile in MinRiverBasins, the region is "camels_mr".
    year
        we use GEE code to generate data for each year, so each year for each region has one data file.
    time_zone
        local time zone and the default is Asia/Hong_Kong (UTC+8)
        Generally our data's time zone is UTC; when we need local time, we need to transform it

    Returns
    -------
    None
    """
    # you can add features or delete features, or change the order, which depends on your txt content
    avg_dataset = [
        "gage_id",
        "time_start",
        "dewpoint_temperature_2m_mean",
        "temperature_2m_mean",
        "skin_temperature_mean",
        "soil_temperature_level_1_mean",
        "soil_temperature_level_2_mean",
        "soil_temperature_level_3_mean",
        "soil_temperature_level_4_mean",
        "lake_bottom_temperature_mean",
        "lake_ice_depth_mean",
        "lake_ice_temperature_mean",
        "lake_mix_layer_depth_mean",
        "lake_shape_factor_mean",
        "lake_total_layer_temperature_mean",
        "snow_albedo_mean",
        "snow_cover_mean",
        "snow_density_mean",
        "snow_depth_mean",
        "snow_depth_water_equivalent_mean",
        "temperature_of_snow_layer_mean",
        "skin_reservoir_content_mean",
        "volumetric_soil_water_layer_1_mean",
        "volumetric_soil_water_layer_2_mean",
        "volumetric_soil_water_layer_3_mean",
        "volumetric_soil_water_layer_4_mean",
        "forecast_albedo_mean",
        "u_component_of_wind_10m_mean",
        "v_component_of_wind_10m_mean",
        "surface_pressure_mean",
        "leaf_area_index_high_vegetation_mean",
        "leaf_area_index_low_vegetation_mean",
    ]
    sum_dataset = [
        "gage_id",
        "time_start",
        "snowfall_hourly_sum",
        "snowmelt_hourly_sum",
        "surface_latent_heat_flux_hourly_sum",
        "surface_net_solar_radiation_hourly_sum",
        "surface_net_thermal_radiation_hourly_sum",
        "surface_sensible_heat_flux_hourly_sum",
        "surface_solar_radiation_downwards_hourly_sum",
        "surface_thermal_radiation_downwards_hourly_sum",
        "evaporation_from_bare_soil_hourly_sum",
        "evaporation_from_open_water_surfaces_excluding_oceans_hourly_sum",
        "evaporation_from_the_top_of_canopy_hourly_sum",
        "evaporation_from_vegetation_transpiration_hourly_sum",
        "potential_evaporation_hourly_sum",
        "runoff_hourly_sum",
        "snow_evaporation_hourly_sum",
        "sub_surface_runoff_hourly_sum",
        "surface_runoff_hourly_sum",
        "total_evaporation_hourly_sum",
        "total_precipitation_hourly_sum",
    ]
    camels_format_index = [
        "Year",
        "Mnth",
        "Day",
        "Hr",
        "dewpoint_temperature_2m",
        "temperature_2m",
        "skin_temperature",
        "soil_temperature_level_1",
        "soil_temperature_level_2",
        "soil_temperature_level_3",
        "soil_temperature_level_4",
        "lake_bottom_temperature",
        "lake_ice_depth",
        "lake_ice_temperature",
        "lake_mix_layer_depth",
        "lake_shape_factor",
        "lake_total_layer_temperature",
        "snow_albedo",
        "snow_cover",
        "snow_density",
        "snow_depth",
        "snow_depth_water_equivalent",
        "temperature_of_snow_layer",
        "skin_reservoir_content",
        "volumetric_soil_water_layer_1",
        "volumetric_soil_water_layer_2",
        "volumetric_soil_water_layer_3",
        "volumetric_soil_water_layer_4",
        "forecast_albedo",
        "u_component_of_wind_10m",
        "v_component_of_wind_10m",
        "surface_pressure",
        "leaf_area_index_high_vegetation",
        "leaf_area_index_low_vegetation",
        "snowfall",
        "snowmelt",
        "surface_latent_heat_flux",
        "surface_net_solar_radiation",
        "surface_net_thermal_radiation",
        "surface_sensible_heat_flux",
        "surface_solar_radiation_downwards",
        "surface_thermal_radiation_downwards",
        "evaporation_from_bare_soil",
        "evaporation_from_open_water_surfaces_excluding_oceans",
        "evaporation_from_the_top_of_canopy",
        "evaporation_from_vegetation_transpiration",
        "potential_evaporation",
        "runoff",
        "snow_evaporation",
        "sub_surface_runoff",
        "surface_runoff",
        "total_evaporation",
        "total_precipitation",
    ]

    gage_id_key = "gage_id"

    # because this function only work for one year and one region, it's better to chose avg and sum files at first
    for f_name in os.listdir(era5_land_dir):
        if fnmatch.fnmatch(
            f_name, "era5_land_" + region + "_avg_mean_" + str(year) + ".csv"
        ):
            avg_data_file = os.path.join(era5_land_dir, f_name)
        if fnmatch.fnmatch(
            f_name, "era5_land_" + region + "_sum_mean_" + str(year) + ".csv"
        ):
            sum_data_file = os.path.join(era5_land_dir, f_name)
    avg_data_temp = pd.read_csv(avg_data_file, sep=",", dtype={avg_dataset[0]: str})
    sum_data_temp = pd.read_csv(sum_data_file, sep=",", dtype={sum_dataset[0]: str})
    # transform UTC time to local
    # refer to: https://stackoverflow.com/questions/33204763/how-to-pass-multiple-arguments-to-the-apply-function
    avg_data_temp[avg_dataset[1]] = avg_data_temp[avg_dataset[1]].apply(
        lambda dt_tmp: utc_to_local(dt_tmp, local_tz=time_zone)
    )
    sum_data_temp[sum_dataset[1]] = sum_data_temp[sum_dataset[1]].apply(
        lambda dt_tmp: utc_to_local(dt_tmp, local_tz=time_zone)
    )
    for i_basin in range(len(gage_dict[gage_id_key])):
        # The next two pieces of code are matching the given dictionary and saving it.It needs to be written in advance
        # The operation is to take the inner txt file corresponding to gage_dict[gage_id_key][i_basin]
        # and save the gage_id in it into the id number of the two input files avg and sum
        # After saving, you can directly run the code, gage_id and name stored in the shp file, you can use GIS to open its properties to view
        # note that the order of gage_id must be aligned with the csv file
        avg_basin_data = avg_data_temp[
            avg_data_temp[avg_dataset[0]].values.astype(str)
            == str(gage_dict[gage_id_key][i_basin])
        ]
        sum_basin_data = sum_data_temp[
            sum_data_temp[avg_dataset[0]].values.astype(str)#np.int64
            == str(gage_dict[gage_id_key][i_basin])
        ]

        if avg_basin_data.shape[0] == 0 or sum_basin_data.shape[0] == 0:
            raise ArithmeticError("chosen basins' number is zero")
        # get Year,Month,Day,Hour info
        csv_date = pd.to_datetime(avg_basin_data[avg_dataset[1]])
        # the hour is set to 12, as 12 is the average hour of a day
        year_month_day_hour = pd.DataFrame(
            [[dt.year, dt.month, dt.day, 12] for dt in csv_date],
            columns=camels_format_index[0:4],
        )
        avg_data_df = pd.DataFrame(
            avg_basin_data.iloc[:, 2:].values,
            columns=camels_format_index[4 : 4 + len(avg_dataset) - 2],
        )
        sum_data_df = pd.DataFrame(
            sum_basin_data.iloc[:, 2:].values,
            columns=camels_format_index[-(len(sum_dataset) - 2) :],
        )
        # concat
        new_data_df = pd.concat([year_month_day_hour, avg_data_df, sum_data_df], axis=1)
        # output the result
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)
        output_file = os.path.join(
            output_dir, gage_dict[gage_id_key][i_basin] + "_lump_era5_land_forcing.txt"
        )
        print(
            "output forcing data of", gage_dict[gage_id_key][i_basin], "year", str(year)
        )
        if os.path.isfile(output_file):
            data_old = pd.read_csv(output_file, sep=" ")
            years = np.unique(data_old[camels_format_index[0]].values)
            if year in years:
                continue
            else:
                os.remove(output_file)
                new_data_df = pd.concat([data_old, new_data_df]).sort_values(
                    by=camels_format_index[0:3]
                )
        new_data_df.to_csv(output_file, header=True, index=False, sep=" ")


