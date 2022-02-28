import pandas as pd
from utils import *

"""
Based on MODIS IGBP, calculate root distribution for a basin (Zeng 2001)

Reference:
Zeng, X. (2001). "Global vegetation root distribution for land modeling." Journal of Hydrometeorology 2(5): 525-530.
https://lpdaac.usgs.gov/products/mcd12q1v006/

Requirement:
(1) processed_igbp.tif: converted IGBP classification in raster form
(2) calculated_root_depth.txt: calculated root_depth 50/99 for each type of land cover based on Eq. (2) and Table 2 in (Zeng 2001)
(2) Catchment shapefiles
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ...
"""


def modis_land_cover_igbp_number2name(index: int):
    names = [
        "EvergreenNeedleleafTree",
        "EvergreenBroadleafTree",
        "DeciduousNeedleleafTree",
        "DeciduousBroadleafTree",
        "MixedForest",
        "ClosedShrubland",
        "OpenShrubland",
        "WoodySavanna",
        "Savanna",
        "Grassland",
        "PermanentWetland",
        "Cropland",
        "UrbanAndBuiltupLand",
        "CroplandOrNaturalVegetaion",
        "SnowAndIce",
        "Barren",
        "WaterBodies",
    ]
    return names[index - 1]


def modis_land_cover_igbp_name2number(name: str):
    names = [
        "EvergreenNeedleleafTree",
        "EvergreenBroadleafTree",
        "DeciduousNeedleleafTree",
        "DeciduousBroadleafTree",
        "MixedForest",
        "ClosedShrubland",
        "OpenShrubland",
        "WoodySavanna",
        "Savanna",
        "Grassland",
        "PermanentWetland",
        "Cropland",
        "UrbanAndBuiltupLand",
        "CroplandOrNaturalVegetaion",
        "SnowAndIce",
        "Barren",
        "WaterBodies",
    ]
    return names.index(name)


class DepthMapper:
    def __init__(self, root_depth_file: str):
        self.land_root_depth = pd.read_table(root_depth_file, sep=",")
        # trans to dict for quick calculation:
        # https://stackoverflow.com/questions/18695605/python-pandas-dataframe-to-dictionary
        self.land_root_depth_50_dict = dict(
            zip(self.land_root_depth["land"], self.land_root_depth["50"])
        )
        self.land_root_depth_99_dict = dict(
            zip(self.land_root_depth["land"], self.land_root_depth["99"])
        )

    def igbp2depth50(self, igbp_index: int):
        name = modis_land_cover_igbp_number2name(igbp_index)
        return self.land_root_depth_50_dict[name]

    def igbp2depth99(self, igbp_index: int):
        name = modis_land_cover_igbp_number2name(igbp_index)
        return self.land_root_depth_99_dict[name]


def root_depth_50_99_stats(
    shape_file: str, igbp_tif: str, depth_mapper: DepthMapper, valid_max=20
):
    """
    the arithmetic mean of catchment effective rooting depth for root_fraction_percentiles=50/99
    For a shapefile, according to IGBP classification, caluculate effective root depth (root_fraction_percentiles=50/99)
    of each grid, and make statistics

    Parameters
    ----------
    shape_file
        shapefile for a basin
    igbp_tif
        converted IGBP classification in raster form
    depth_mapper
        DepthMapper object; map IGBP to effective depth
    valid_max
        there are 17 classes, here we set valid_max = 20 just to exclude invalid values;
        invalid value is typically 255 because the data type is uint8 (range: 0 - 255),
        also exists some strange case, for instance, mine is 241

    Returns
    -------
    dict
        {'root_depth_50': np.mean(depth50), 'root_depth_99': np.mean(depth99)}
    """
    res = extract_raster_by_shape_file(
        raster=igbp_tif, shape_file=shape_file, output_file=None
    )
    res = res[res != -9999]
    res_list = res[res < valid_max].flatten().tolist()
    # mapping igbp classification to effective rooting depth for each pixel
    depth50 = [depth_mapper.igbp2depth50(index) for index in res_list]
    depth99 = [depth_mapper.igbp2depth99(index) for index in res_list]
    return {"root_depth_50": np.mean(depth50), "root_depth_99": np.mean(depth99)}
