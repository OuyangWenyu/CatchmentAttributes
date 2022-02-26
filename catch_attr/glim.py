import numpy as np
import pandas as pd

from utils import extract_raster_by_shape_file

"""
According to the GliM dataset, calculate geology attributes.
Reference: Hartmann, J., Moosdorf, N., 2012. The new global lithological map database GLiM: A representation of rock 
properties at the Earth surface. Geochemistry, Geophysics, Geosystems, 13. DOI: 10.1029/2012GC004370

Requirement: 
GlimRaster.tif: converted GliM in raster form (you can use ArcMap/QGis to complete this operation, we cannot directly share the converted data due to copyright issues)
GLiMCateNumberMapping.csv: mapping number to lithology category
glim_short_long_name.txt: mapping short name to long name of lithology categories
"""


class Glim:
    def __init__(
        self,
        glim_raster_tif: str,
        glim_cate_number_mapping_file: str,
        short2long_name_txt: str,
        nan_value=-9999,
    ):
        """
        Calculate basins' lithology ratio based on GliM data set

        Parameters
        ----------
        glim_raster_tif
            Glim tif file
        glim_cate_number_mapping_file
            A csv file: map from GLiM raster LithoValue (raster value) to Litho (litho class);
        short2long_name_txt
            A txt file: map from litho short name to long name; glim_short2longname.txt
        nan_value
            default is -9999
        """
        self.glim_raster_tif = glim_raster_tif
        self.short2long_dataframe = pd.read_table(short2long_name_txt, sep=",")
        self.glim_mapping_dataframe = pd.read_table(
            glim_cate_number_mapping_file, sep=","
        )
        # self.glim_mapping_dataframe['xx'] = [s[:2] for s in self.glim_mapping_dataframe['Litho']]
        self.nan_value = nan_value

    def glim_number2geol_mapping(self, value: int):
        return self.glim_mapping_dataframe[
            self.glim_mapping_dataframe["xxValue"] == value
        ]["xx"].values[0]

    def glim_geol2number_mapping(self, geol: str):
        return self.glim_mapping_dataframe[self.glim_mapping_dataframe["xx"] == geol][
            "xxValue"
        ].values

    def short2long_name(self, short_name: str):
        return self.short2long_dataframe[
            self.short2long_dataframe["short"] == short_name
        ]["long"].values[0]

    def extract_basin_attributes_glim_all(self, shape_file: str) -> dict:
        """
        Extract geology attributes for one basin from GLiM

        Parameters
        ----------
        shape_file
            shapefile path
        Returns
        -------
        dict
            ratio of each geology types
        """
        res = extract_raster_by_shape_file(
            raster=self.glim_raster_tif, shape_file=shape_file, output_file=None
        )
        res = res[res < 1000].flatten()
        res_list = res[res != self.nan_value].flatten().tolist()
        geol_class_number, count = np.unique(res_list, return_counts=True)
        geol_class = [
            self.glim_number2geol_mapping(number) for number in geol_class_number
        ]
        res = {}
        for name, c in zip(geol_class, count):
            res[name] = c / np.sum(count)
        return res

    def extract_basin_attributes_glim(self, shape_file: str) -> dict:
        """
        Extract geology attributes for one basin from GLiM

        Parameters
        ----------
        shape_file
            shapefile path
        Returns
        -------
        dict
            geology attributes
        """
        res = extract_raster_by_shape_file(
            raster=self.glim_raster_tif, shape_file=shape_file, output_file=None
        )
        res = res[res < 1000].flatten()
        res_list = res[res != self.nan_value].flatten().tolist()
        geol_class_number, count = np.unique(res_list, return_counts=True)
        geol_class = np.array(
            [self.glim_number2geol_mapping(number) for number in geol_class_number]
        )
        geol_class_rank = [x for _, x in sorted(zip(count, geol_class), reverse=True)]
        if len(geol_class_rank) == 0:
            return {
                "geol_class_1st": None,
                "geol_class_1st_frac": None,
                "geol_class_2nd": None,
                "geol_class_2nd_frac": None,
                "carb_rocks_frac": None,
            }
        geol_class_1st = geol_class_rank[0]
        geol_class_1st_count = count[geol_class == geol_class_1st]
        geol_class_1st_frac = (geol_class_1st_count / np.sum(count))[0]
        if len(geol_class_rank) > 1:
            geol_class_2nd = geol_class_rank[1]
            geol_class_2nd_count = count[geol_class == geol_class_2nd]
            geol_class_2nd_frac = (geol_class_2nd_count / np.sum(count))[0]
        else:
            geol_class_2nd = None
            geol_class_2nd_frac = 0

        carb_rocks_count = count[geol_class == "sc"]
        if len(carb_rocks_count) == 0:
            carb_rocks_frac = 0
        else:
            carb_rocks_frac = (carb_rocks_count / np.sum(count))[0]

        return {
            "geol_class_1st": geol_class_1st,
            "geol_class_1st_frac": geol_class_1st_frac,
            "geol_class_2nd": geol_class_2nd,
            "geol_class_2nd_frac": geol_class_2nd_frac,
            "carb_rocks_frac": carb_rocks_frac,
        }
