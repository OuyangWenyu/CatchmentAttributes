import os
from qgis import processing

output_dir = os.path.join(
    "D:\\", "code", "CatchmentAttributes", "data", "soil_source_data"
)
input_files = ["CLYPPT_M_sl" + str(i) + "_250m_ll" for i in range(1, 8)]
out_files = [
    os.path.join(output_dir, input_file + "_downscaled.tif")
    for input_file in input_files
]
for i in range(len(out_files)):
    results = processing.runAndLoadResults(
        "gdal:warpreproject",
        {
            "INPUT": input_files[i],
            "RESAMPLING": 5,
            "NODATA": -9999,
            "TARGET_RESOLUTION": 0.05,
            "OUTPUT": out_files[i],
        },
    )
