from utils import *

"""
Basedo n MODIS MCD12Q1 product LC_Type1, calculate fractions of each land cover in a basin
(Annual International Geosphere-Biosphere Programme (IGBP) classification)

Reference:
https://lpdaac.usgs.gov/products/mcd12q1v006/

Requirement:
(1) processed_igbp.tif: converted IGBP classification in raster form MODIS mcd12q1v006
(2) Catchment shapefiles
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ......
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
    try:
        return names[index - 1]
    except IndexError:
        return "nan"


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


def igbp_stats(shapefile: str, igbp_tif: str, valid_max=20):
    """
    For the basin in a shapefile, calculate dom_land_cover, dom_land_cover_frac, forest_frac according to Modis_IGBP

    Parameters
    ----------
    shapefile
        the basin's shapefile
    igbp_tif
        processed_igbp.tif file
    valid_max
        there are 17 classes, here we set valid_max = 20 just to exclude invalid values;
        invalid value is typically 255 because the data type is uint8 (range: 0 - 255),
        also exists some strange case, for instance, mine is 241

    Returns
    -------
    dict
        {'dom_land_cover': land_class_1st, 'dom_land_cover_frac': land_class_1st_frac, 'forest_frac': forest_frac}
    """

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

    res = extract_raster_by_shape_file(
        raster=igbp_tif, shape_file=shapefile, output_file=None
    )
    res = res[res != -9999].flatten()
    res_list = res[res < valid_max].flatten().tolist()
    res_str = [modis_land_cover_igbp_number2name(number) for number in res_list]
    land_class, count = np.unique(res_str, return_counts=True)
    # land_class_rank = [x for _, x in sorted(zip(count, land_class), reverse=True) if x != 'nan']

    res = {}
    for name in names:
        res[name] = 0
    for num, name in zip(count, land_class):
        if name != "nan":
            res[name] = num / np.sum(count)

    print("\nshapefile:", shapefile)
    print(res)
    return res
