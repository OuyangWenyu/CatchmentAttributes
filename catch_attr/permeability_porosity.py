from utils import *

"""
Calculate catchment aggregated permeability and porosity based on GLHYMPS.

Data:
Source: https://dataverse.scholarsportal.info/dataset.xhtml?persistentId=doi:10.5683/SP2/DLGXYO
Processed: https://zenodo.org/record/5137288/files/9_code_data.zip?download=1

Reference: 
Gleeson, Tom, 2018, "GLobal HYdrogeology MaPS (GLHYMPS) of permeability and porosity", https://doi.org/10.5683/SP2/DLGXYO, Scholars Portal Dataverse, V1

"""


class GLHYMPS:
    def __init__(
        self,
        permeabilit_no_permafrost_raster_tif: str,
        porosity_raster_tif: str,
        nan_value=65535,
    ):
        self.permeabilit_no_permafrost_raster_tif = permeabilit_no_permafrost_raster_tif
        self.porosity_raster_tif = porosity_raster_tif
        self.nan_value = nan_value

    def zonal_stats_glhymps(self, shape_file: str) -> dict:
        permeability = zonal_stats_singletif(
            tif_file=self.permeabilit_no_permafrost_raster_tif,
            shape_file=shape_file,
            valid_min=-self.nan_value,
            valid_max=self.nan_value,
        )
        porosity = zonal_stats_singletif(
            tif_file=self.porosity_raster_tif,
            shape_file=shape_file,
            valid_min=-self.nan_value,
            valid_max=self.nan_value,
        )

        return {"permeability": permeability, "porosity": porosity}


def shp_id(shpfile: str):
    return re.findall(r"[\d]+", shpfile)[-1]
