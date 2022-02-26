import os
import numpy as np
import pandas as pd
import pytest

from catch_attr.climate import (
    p_mean,
    high_prec_freq,
    high_prec_dur,
    high_prec_timing,
    low_prec_freq,
    low_prec_dur,
    low_prec_timing,
    frac_snow_daily,
)

import definitions
from catch_attr.basin_era5_process import utc_to_local, trans_era5_land_to_camels_format


@pytest.fixture
def forcing_data():
    forcing_file = os.path.join(
        definitions.DATASET_DIR, "61019_lump_era5_land_forcing.txt"
    )
    df = pd.read_csv(forcing_file, sep="\s+")
    return df


def test_tz_trans():
    utc_time_str = "2020-01-01T00:00:00"
    local_tz = "Asia/Hong_Kong"
    utc_ts = utc_to_local(utc_time_str, local_tz)
    assert utc_ts == "2020-01-01T08:00:00"


def test_tz_np_trans():
    utc_time = np.datetime64("2020-01-01T00:00:00")
    utc_ts = utc_to_local(utc_time)
    assert utc_ts == "2020-01-01T08:00:00"


def test_gee_daily_era5_land_to_camels_format():
    output_dir = os.path.join("results", "basin_mean_forcing")
    region = "camels_mr"
    year = 2010
    gage_dict = pd.read_csv(
        os.path.join(definitions.DATASET_DIR, "camels_mr_name.txt"),
        dtype={"gage_id": str},
    ).to_dict(orient="list")
    trans_era5_land_to_camels_format(
        definitions.DATASET_DIR, output_dir, gage_dict, region, year
    )
    print("Trans finished")


def test_p_eman(forcing_data):
    p_mean_value = p_mean(forcing_data["total_precipitation"])
    assert p_mean_value == 0.003689500707409431


def test_high_prec_freq(forcing_data):
    high_prec_freq_ = high_prec_freq(forcing_data["total_precipitation"])
    assert high_prec_freq_ == 4.425107547907705


def test_high_prec_dur(forcing_data):
    high_prec_dur_ = high_prec_dur(forcing_data["total_precipitation"])
    assert high_prec_dur_ == 1.1923076923076923


def test_high_prec_timing(forcing_data):
    high_prec_timing_ = high_prec_timing(forcing_data)
    assert high_prec_timing_ == "jja"


def test_low_prec_freq(forcing_data):
    low_prec_freq_ = low_prec_freq(forcing_data["total_precipitation"])
    assert low_prec_freq_ == 116.05201407899882


def test_low_prec_dur(forcing_data):
    low_prec_dur_ = low_prec_dur(forcing_data["total_precipitation"])
    assert low_prec_dur_ == 2.862676056338028


def test_low_prec_timing(forcing_data):
    low_prec_timing_ = low_prec_timing(forcing_data)
    assert low_prec_timing_ == "djf"


def test_frac_snow_daily(forcing_data):
    frac_snow_daily_ = frac_snow_daily(forcing_data)
    assert frac_snow_daily_ == 0.4106374657802112
