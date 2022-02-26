import math
import richdem as rd
import shapefile

from utils import *

"""
According to ASTER GDEM: https://asterweb.jpl.nasa.gov/gdem.asp , make statistics for topography of a basin

Requirement:
(1) ASTER GDEM
├── DEM
|   ├── ASTGTMV003_N34E111_dem.tif
|   ├── ASTGTMV003_N32E110_dem.tif
|   ├── ASTGTMV003_N33E109_dem.tif
|   ├── ASTGTMV003_N34E108_dem.tif
|   ├── ...
(2) Catchment shapefiles
├── shapefiles
|   ├── basin_xxxx.shp
|   ├── ...
"""


def load_n_e_from_dem_name(dem_file: str):
    """
    Get lat and lon from aster dem file names

    Parameters
    ----------
    dem_file


    Returns
    -------

    """
    n = re.findall("N(\d+)", dem_file)[0]
    e = re.findall("E(\d+)", dem_file)[0]
    return {"N": int(n), "E": int(e)}


def shapefile_n_e(shpfile: str):
    """
    Get min/max lat/lon, this is for determining the range of needed dem files

    Parameters
    ----------
    shpfile
        a basin's shapefile

    Returns
    -------
    dict
        box of the basin polygon -- {'N_min': bbox[1], 'N_max': bbox[3], 'E_min': bbox[0], 'E_max': bbox[2]}
    """
    sf = shapefile.Reader(shpfile)
    bbox = sf.bbox
    return {"N_min": bbox[1], "N_max": bbox[3], "E_min": bbox[0], "E_max": bbox[2]}


def fetch_shapefile_needed_dem_range(shpfile: str):
    """
    Get the range of needed dem files for the given shapefile

    Parameters
    ----------
    shpfile
        a basin's shapefile

    Returns
    -------
    dict
        dem's range -- {'Ns': ns, 'Es': es}
    """
    n_e = shapefile_n_e(shpfile)
    ns = range(math.floor(n_e["N_min"]), math.ceil(n_e["N_max"]) + 1)
    es = range(math.floor(n_e["E_min"]), math.ceil(n_e["E_max"]) + 1)
    return {"Ns": ns, "Es": es}


def merge_and_reproject_dems(
    shpfile: str, dem_folder: str, tmp_merged: str, tmp_reprojected: str
):
    """
    Select required DEM files, merge them and reproject the merged tif file; then, save it to "TOPO" directory

    Parameters
    ----------
    shpfile
        a basin's shapefile
    dem_folder:
        directory where we put DEM tiles
    tmp_merged
        merged DEM tif file path
    tmp_reprojected
        reprojected merged DEM tif file path

    Returns
    -------
    None
    """
    if os.path.isfile(tmp_merged):
        os.remove(tmp_merged)
    if os.path.isfile(tmp_reprojected):
        os.remove(tmp_reprojected)
    # get needed N E range
    needed_ns = fetch_shapefile_needed_dem_range(shpfile)["Ns"]
    needed_es = fetch_shapefile_needed_dem_range(shpfile)["Es"]

    # get needed tifs
    files = absolute_file_paths(dem_folder)
    files = [file for file in files if file.endswith("dem.tif")]
    needed_tifs = []

    for file in files:
        n = load_n_e_from_dem_name(file)["N"]
        e = load_n_e_from_dem_name(file)["E"]
        if (n in needed_ns) and (e in needed_es):
            needed_tifs.append(file)

    # merge tifs
    if len(needed_tifs) == 0:
        raise FileNotFoundError(
            f"did not find needed tifs for determining topograpy attributes | shpfile: {shpfile}"
        )
    merge_tifs(needed_tifs, tmp_merged)
    reproject_tif(tmp_merged, tmp_reprojected)


def elev_mean(shpfile: str, tmp_reprojected: str) -> float:
    """
    Calculate mean elevation of the catchment

    Parameters
    ----------
    shpfile
        a basin's shapefile
    tmp_reprojected
        reprojected merged DEM tif file path

    Returns
    -------
    float
        mean elevation of the basin
    """
    # zonal stats
    try:
        res = zonal_stats_singletif(tmp_reprojected, shpfile)
    except ValueError as e:
        print(e)
        print("Please merge required DEMs, reproject the merged one and save it!!!")
        return np.nan
    return res


def calculate_slope(dem):
    """
    Calculate the slope of a given dem

    Parameters
    ----------
    dem
        dem file path

    Returns
    -------
    np.array
        slope
    """
    dem_path = dem
    shasta_dem = rd.LoadGDAL(dem_path, no_data=-128)
    slope = rd.TerrainAttribute(shasta_dem, attrib="slope_riserun")
    return np.array(slope.data)


def slope_mean(shpfile: str, tmp_reprojected: str, tmp_slope: str):
    """
    Calculate the slope of a given catchment

    Parameters
    ----------
    shpfile
        a basin's shapefile
    tmp_reprojected
        reprojected merged DEM tif file path
    tmp_slope
        slope file

    Returns
    -------
    float
        mean slope of a basin
    """
    # get slope
    slope = calculate_slope(tmp_reprojected)

    # rise (m, dem) / run (degree, coordinates)
    slope = slope / 111  # 1 degree = 111km

    # get raster range
    ds = gdal.Open(tmp_reprojected)
    width = ds.RasterXSize
    height = ds.RasterYSize
    gt = ds.GetGeoTransform()
    minx = gt[0]
    miny = gt[3] + width * gt[4] + height * gt[5]
    maxx = gt[0] + width * gt[1] + height * gt[2]
    maxy = gt[3]
    ds = None

    # write slope tif
    geotif_from_array(
        slope,
        lat_start=miny,
        lat_end=maxy,
        lon_start=minx,
        lon_end=maxx,
        degree=0.1,
        output_file=tmp_slope,
    )

    # zonal stats
    try:
        res = zonal_stats_singletif(tmp_slope, shpfile)
    except ValueError as e:
        print(e)
        print("Please merge required DEMs, reproject the merged one and save it!!!")
        return np.nan
    return res
