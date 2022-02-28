import shapefile
from tqdm import tqdm
from shapely.geometry import Point, LineString, Polygon
from functools import partial
import shapely.ops as ops
import pyproj
from utils import *
from shapely.ops import transform

"""
Drainage basin boundary data and the river network data are obtained from the Global Drainage Basin Database (GDBD) 
dataset: https://www.cger.nies.go.jp/db/gdbd/gdbd_index_e.html. Here, determining the basin outlet needs river network 
and basin boundaries as input. Since the river network provided by GDBD did not cover all basins (mainly watersheds 
where the river stream is not clear), basin outlet and hence river length cannot be derived for some of the basins.
"""


def latlon2km(p1, p2):
    """
    Convert distance in the degree to km

    Parameters
    ----------
    p1
        point 1
    p2
        point 2

    Returns
    -------
    float
        distance (km)
    """
    pyproj.Proj("+init=epsg:4326")
    line = LineString([p1, p2])
    wgs84 = pyproj.Proj(init="epsg:4326")
    utm = pyproj.Proj(init="epsg:32649")

    project = partial(pyproj.transform, wgs84, utm)

    utm_polyline = transform(project, line)
    length_km = utm_polyline.length / 1000
    return length_km


def find_outlet(catchment_shp, stream_shps):
    """
    Find catchment outlet point given river stream shps and catchment shapefile

    Parameters
    ----------
    catchment_shp
        shapefile of a basin
    stream_shps
        shapefile of streams

    Returns
    -------
    tuple
        coordination of outlet
    """
    basin = shapefile.Reader(catchment_shp).shapeRecord(0).shape
    stream_shapes = shapefile.Reader(stream_shps).shapes()

    basin_polygon = Polygon(basin.points)
    if not basin_polygon.is_valid:
        print("invalid polygon")
        return
    for stream in stream_shapes:
        line = LineString(stream.points)
        intersection = basin_polygon.exterior.intersection(line)
        if intersection.is_empty:
            continue
        else:
            try:
                return intersection.coords.xy[0][0], intersection.coords.xy[1][0]
            except:
                continue


def get_record(basins, i):
    """
    Read shapefile information

    Parameters
    ----------
    basins
        shapefile.Reader(shapefile)
    i
        return i-th polygon's information

    Returns
    -------
    dict
        ith basin's information {"key1":value1,"key2":value2,...}
    """
    fields = basins.fields[1:]
    field_names = [field[0] for field in fields]
    atr = dict(zip(field_names, basins.shapeRecord(i).record))
    return atr


def catchment_perimeter(catchment_shp):
    """
    Calculate the Perimeter of the catchment given shapefile

    Parameters
    ----------
    catchment_shp
        shapefile of a basin

    Returns
    -------
    float
        perimeter of the basin (km)
    """
    polyline = Polygon(shapefile.Reader(catchment_shp).shapeRecord(0).shape.points)
    wgs84 = pyproj.Proj(init="epsg:4326")
    utm = pyproj.Proj(init="epsg:32649")

    project = partial(pyproj.transform, wgs84, utm)

    utm_polyline = transform(project, polyline)
    length_km = utm_polyline.length / 1000
    return length_km


def longest_distance(catchment_shp, stream_shps):
    """
    For a given point, find the remotest point on the polygon boundary and return the distance

    Parameters
    ----------
    catchment_shp
        shapefile of a basin
    stream_shps
        shapefile of streams

    Returns
    -------
    float
        longest distance of a basin to its outlet
    """
    basin = shapefile.Reader(catchment_shp).shapeRecord(0).shape
    outlet = find_outlet(catchment_shp, stream_shps)
    if outlet:
        max_dis = 0
        for p in basin.points:
            if Point(p).distance(Point(outlet)) > max_dis:
                max_dis = Point(p).distance(Point(outlet))
                max_point = p
        return latlon2km(max_point, outlet)


def form_factor(area, length):
    """
    Catchment Form factor

    [Section Catchment Characteristics] Subramanya, K. (2013). Engineering hydrology, 4e. Tata McGraw-Hill Education.

    Parameters
    ----------
    area
        area of a basin
    length
        longest distance of a basin to its outlet

    Returns
    -------
    float
        catchment form factor
    """
    return area / length**2


def shape_factor(area, length):
    """
    Catchment Shape factor

    [Section Catchment Characteristics] Subramanya, K. (2013). Engineering hydrology, 4e. Tata McGraw-Hill Education.

    Parameters
    ----------
    area
        area of a basin
    length
        longest distance of a basin to its outlet

    Returns
    -------
    float
        catchment shape factor
    """

    return length**2 / area


def compactness_coefficient(perimeter, area):
    """
    Catchment Compactness coefficient

    [Section Catchment Characteristics] Subramanya, K. (2013). Engineering hydrology, 4e. Tata McGraw-Hill Education.

    Parameters
    ----------
    perimeter
        perimeter of a basin
    area
        area of a basin

    Returns
    -------
    float
        catchment compactness coefficient
    """
    return 0.2821 * perimeter / np.sqrt(area)


def circulatory_ratio(perimeter, area):
    """
    Catchment Circulatory ratio

    [Section Catchment Characteristics] Subramanya, K. (2013). Engineering hydrology, 4e. Tata McGraw-Hill Education.

    Parameters
    ----------
    perimeter
        perimeter of a basin
    area
        area of a basin

    Returns
    -------
    float
        catchment circulatory ratio
    """
    return 12.57 * area / perimeter**2


def elongation_ratio(area, length):
    """
    Catchment Elongation ratio

    [Section Catchment Characteristics] Subramanya, K. (2013). Engineering hydrology, 4e. Tata McGraw-Hill Education.

    Parameters
    ----------
    area
        area of a basin
    length
        longest distance of a basin to its outlet

    Returns
    -------
    float
        catchment elongation ratio
    """
    return 1.128 * np.sqrt(area) / length


def basin_area(basin_shp):
    """
    Calculate catchment area given shapefile

    Parameters
    ----------
    basin_shp
        shapefile of a basin

    Returns
    -------
    float
        area of the basin (km^2)
    """
    basin = shapefile.Reader(basin_shp).shapeRecord(0).shape
    geom = Polygon(basin.points)
    geom_area = ops.transform(
        partial(
            pyproj.transform,
            pyproj.Proj(init="EPSG:4326"),
            pyproj.Proj(proj="aea", lat_1=geom.bounds[1], lat_2=geom.bounds[3]),
        ),
        geom,
    )
    return geom_area.area / 1000**2  # km^2


def basin_topo_stats(basin_shps, stream_shps):
    """
    Catchment shape factors statistics given a list of basin shapefiles and river stream shapefiles

    Parameters
    ----------
    basin_shps
        shapefiles of all target basins
    stream_shps
        shapefile of streams

    Returns
    -------
    dict
        catchment shape factors statistics
    """
    res = {}
    for basin_shp in tqdm(basin_shps):
        basin = shapefile.Reader(basin_shp)
        gage_id = int(get_record(basin, 0)["ID"])
        ld = longest_distance(basin_shp, stream_shps)
        ba = basin_area(basin_shp)
        cp = catchment_perimeter(basin_shp)
        if ld:
            res[gage_id] = {
                "Length": ld,
                "Area": ba,
                "FormFactor": form_factor(ba, ld),
                "ShapeFactor": shape_factor(ba, ld),
                "CompactnessCoefficient": compactness_coefficient(cp, ba),
                "CirculatoryRatio": circulatory_ratio(cp, ba),
                "ElongationRatio": elongation_ratio(ba, ld),
            }
        # If the given polygon is invalid, the length attribute (L) cannot be determined, and other variables depend
        # on that.
        else:
            res[gage_id] = {
                "Length": None,
                "Area": ba,
                "FormFactor": None,
                "ShapeFactor": None,
                "CompactnessCoefficient": None,
                "CirculatoryRatio": None,
                "ElongationRatio": None,
            }
    return res
