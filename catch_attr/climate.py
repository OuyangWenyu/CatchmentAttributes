from itertools import groupby
import pandas as pd
import numpy as np
from sklearn.linear_model import LinearRegression
import datetime


def p_seasonality(data: pd.DataFrame):
    """
    seasonality and timing of precipitation (estimated using sine curves to represent the annual temperature and
    precipitation cycles; positive (negative) values indicate that precipitation peaks in summer (winter); values close
    to 0 indicate uniform precipitation throughout the year)
    Reference:
    Woods R A. Analytical model of seasonal climate impacts on snow hydrology: Continuous snowpacks[J].
    Advances in Water Resources, 2009, 32(10): 1465-1481.
    Method detail is described at the end of section 2.3 in the original paper.
    The parameters were estimated by exhaustive search on st, combined with least squares regression for tbar and delta_t;
    the same method was used to estimate pbar; delta_p and sp.

    Parameters
    ----------
    data pd.DataFrame containing PRE ['20-20时累计降水量'] and TEM ['平均气温'] columns | longer is better

    Returns
    -------
    float
        p_seasonality
    """
    df_date = data[["Year", "Mnth", "Day"]]
    df_date.columns = ["year", "month", "day"]
    date = pd.to_datetime(df_date).values.astype("datetime64[D]")
    data.index = date
    data = data.loc[datetime.datetime(2019, 1, 1): datetime.datetime(2019, 12, 31)]

    best_score = 0
    for st in range(365):
        X = []
        y = []
        for i in range(len(data)):
            t = data.index[i].day_of_year
            X.append(np.sin(2 * np.pi * (t - st) / 365))
            y.append(data.iloc[i]["temperature_2m"] - 273.15)
        #             y.append(data.iloc[i]['平均气温'])
        X = np.array(X).reshape(-1, 1)
        y = np.array(y)
        reg = LinearRegression(fit_intercept=True, n_jobs=-1).fit(X, y)
        if reg.score(X, y) > best_score:
            best_reg = reg
            best_score = reg.score(X, y)
            best_X = X
            best_st = st
    delta_t = best_reg.coef_[0]
    st = best_st
    best_reg_t = best_reg

    best_score = 0
    for sp in range(365):
        X = []
        y = []
        for i in range(len(data)):
            t = data.index[i].day_of_year
            X.append(np.sin(2 * np.pi * (t - sp) / 365))
            y.append(data.iloc[i]["total_precipitation"] * 1000)
        #             y.append(data.iloc[i]['20-20时累计降水量'])
        X = np.array(X).reshape(-1, 1)
        y = np.array(y)
        reg = LinearRegression(fit_intercept=True, n_jobs=-1).fit(X, y)
        if reg.score(X, y) > best_score:
            best_reg = reg
            best_score = reg.score(X, y)
            best_X = X
            best_sp = sp
    pbar = best_reg.intercept_
    delta_p = best_reg.coef_[0] / pbar
    sp = best_sp
    best_reg_p = best_reg

    p_season = delta_p * np.sign(delta_t) * np.cos(2 * np.pi * (sp - st) / 365)
    return p_season, delta_t, st, delta_p, sp, best_reg_t, best_reg_p


def split_a_list_at_zeros(L):
    return [list(g) for k, g in groupby(L, key=lambda x: x != 0) if k]


def series_mean(data: pd.Series) -> float:
    """
    Mean value of daily variable of a basin

    Parameters
    ----------
    data
        data of one var of a basin

    Returns
    -------
    float
        mean value of the var of a basin
    """
    return float(data.mean())


def high_prec_freq(data: pd.Series) -> float:
    """
    Number of days in which over-5-times-mean-value precipitation happened of one year

    Parameters
    ----------
    data
        precipitation data of a basin

    Returns
    -------
    float
        mean high precipitation days in one year
    """
    num_high_pre_days = len(data[data > data.mean() * 5].dropna())
    return num_high_pre_days / len(data) * 365


def high_prec_dur(data: pd.Series) -> float:
    """
    Mean duration of all high precipitations

    Parameters
    ----------
    data
        precipitation data of a basin

    Returns
    -------
    float
        mean duration of all high precipitations
    """
    data = np.array(data)
    tmp_data = data.copy()
    tmp_data[tmp_data < data.mean() * 5] = 0
    tmp = [len(x) for x in split_a_list_at_zeros(tmp_data)]
    if len(tmp) > 0:
        return np.array(tmp).mean()
    else:
        return 0


def high_prec_timing(data: pd.DataFrame) -> str:
    """
    The season with the most frequent high precipitations

    Parameters
    ----------
    data
        forcing data

    Returns
    -------
    float
        season with the most frequent high precipitations
    """
    months = [
        x
        for x in data[
            data["total_precipitation"] > data["total_precipitation"].mean() * 5
            ]["Mnth"]
    ]
    seasons = [month2season(x) for x in months]
    seasons, counts = np.unique(seasons, return_counts=True)
    if len(counts) > 0:
        return [x for _, x in sorted(zip(counts, seasons))][-1]
    else:
        return ""


def month2season(month: int) -> str:
    """
    Which season of the month

    DJF=Dec-Feb, MAM=Mar-May, JJA=Jun-Aug, SON=Sep-Nov

    Parameters
    ----------
    month
        1-12

    Returns
    -------
    str
        season
    """
    if month in [3, 4, 5]:
        return "mam"
    elif month in [6, 7, 8]:
        return "jja"
    elif month in [9, 10, 11]:
        return "son"
    elif month in [12, 1, 2]:
        return "djf"


def low_prec_freq(data: pd.Series) -> float:
    """
    Number of days in which precipitaion < 1e-3 m/d in one year

    Parameters
    ----------
    data
        precipitation data of a basin

    Returns
    -------
    float
        mean number of low precipitation days in one year
    """
    num_low_pre_days = len(data[data < 1e-3].dropna())
    return num_low_pre_days / len(data) * 365


def low_prec_dur(data: pd.Series) -> float:
    """
    Mean duration of all low precipitations

    Parameters
    ----------
    data
        precipitation data of a basin

    Returns
    -------
    float
        mean duration of all low precipitations
    """
    data = np.array(data)
    tmp_data = data.copy()
    tmp_data[data < 1e-3] = 1
    tmp_data[data > 1e-3] = 0
    tmp = [len(x) for x in split_a_list_at_zeros(tmp_data)]
    if len(tmp) > 0:
        return np.array(tmp).mean()
    else:
        return 0


def low_prec_timing(data: pd.DataFrame) -> str:
    """
    The season with the most frequent high precipitations

    Parameters
    ----------
    data
        forcing data

    Returns
    -------
    float
        season with the most frequent high precipitations
    """
    months = [x for x in data[data["total_precipitation"] < 1e-3]["Mnth"]]
    seasons = [month2season(x) for x in months]
    seasons, counts = np.unique(seasons, return_counts=True)
    return [x for _, x in sorted(zip(counts, seasons))][-1]


def frac_snow_daily(df: pd.DataFrame) -> float:
    """
    The fraction of precipitation falling as snow (for days colder than 0 °C)

    Parameters
    ----------
    df
        forcing data

    Returns
    -------
    float
        fraction of snow days in all days
    """
    snow_days = df[(df["temperature_2m"] < 273.15) & (df["total_precipitation"] > 0)]
    return len(snow_days) / len(df)
