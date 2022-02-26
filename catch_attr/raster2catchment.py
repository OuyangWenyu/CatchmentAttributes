import numpy as np
import random
from datetime import datetime
from multiprocessing import Process
import pandas as pd
from tqdm import tqdm
import os

"""
将插值好的气象栅格数据(使用raster_surf.py)转换为流域的面均值，使用采样法计算。

Requirement:
(1) 插值好的栅格
├── folder_raster
|   ├── 20-20时累计降水量
|   |   ├── 1954-1-1-20-20时累计降水量
|   |   ├── 1954-1-2-20-20时累计降水量
|   ├── 大型蒸发量
|   |   ├── ...
(2) Catchment shapefiles
├── folder_shp
|   ├── outwtrshd_0000.shp
|   ├── outwtrshd_0000.dbf
|   ├── outwtrshd_0000.sbx
|   ├── outwtrshd_0000.cpg
|   ├── ...
"""

# 指定插值范围
lat_start = 15
lat_end = 55
lon_start = 70
lon_end = 140
degree = 0.1

xi = np.round(np.arange(lat_start, lat_end, degree), 1)
yi = np.round(np.arange(lon_start, lon_end, degree), 1)


def load_list(path):
    """

    :param path: .txt 文件路径
    :return: .txt 文件内容
    """
    score = []
    with open(path, "r") as f:
        for line in f:
            score.append(line.strip())
    return score


def load_json(path):
    """

    :param path: .json 文件路径
    :return: .json 文件内容
    """
    import json

    with open(path, "r", encoding="utf8") as read_file:
        data = json.load(read_file)
    return data


def absoluteFilePaths(directory):
    """

    :param directory: 文件夹路径
    :return: 文件夹以及子文件夹所有文件路径
    """
    import os

    def absoluteFilePaths(directory):
        for dirpath, _, filenames in os.walk(directory):
            for f in filenames:
                yield os.path.abspath(os.path.join(dirpath, f))

    return list(absoluteFilePaths(directory))


def read_tif(tif_file):
    """

    :param tif_file: .tif 文件路径
    :return: tif 文件转换成 numpy.array
    """
    from PIL import Image

    im = Image.open(tif_file)
    return np.array(im)


def shp_points(shp):
    """

    :param shp: shp 文件路径
    :return: 将 shp 文件转换为坐标列表
    """
    import shapefile

    return shapefile.Reader(shp).shape(0).points


def tif_shp_index_mean(tif, points, num_sample):
    """

    :param tif: tif 文件路径
    :param points: x,y 点坐标列表
    :param num_sample: 要采样的点的个数，个数越多，采样密度越大
    :return: 计算流域面均
    """
    arr = read_tif(tif)
    if len(points) > num_sample:
        points = random.sample(points, num_sample)
    values = []
    for sample in points:
        y, x = sample
        x_index = np.where(xi == x)
        y_index = np.where(yi == y)
        values.append(arr[x_index, y_index])
    return np.mean(values)


def one_shp(name, num_sample, tifs, shp_points_d, outdir):
    """

    :param name: 流域名称
    :param num_sample: 计算面均采样个数，默认100000（全部采样）
    :pram tifs: tif 文件路径列表（folder_raster 文件夹中包含的 tif 文件列表）
    :param shp_points_d: dict: {name: {[x1, y1], [x2, y2}}
    :param outdir: 输出路径
    :return: 统计给定 shapefile 在所有气象插值栅格（tifs）上的面均值
    """
    res = {}
    for tif in tifs:
        if "降水量" in tif:
            year, month, day = tuple(tif.split("\\")[-1].split(".")[0].split("-"))[:3]
            var = "20-20时累计降水量"
        else:
            year, month, day, var = tuple(tif.split("\\")[-1].split(".")[0].split("-"))
        year = int(year)
        month = int(month)
        day = int(day)

        if datetime(year, month, day) not in res:
            res[datetime(year, month, day)] = {}
        res[datetime(year, month, day)][var] = tif_shp_index_mean(
            tif, shp_points_d[name], num_sample=num_sample
        )
    if not os.path.isdir(f"{outdir}/{name}"):
        os.mkdir(f"{outdir}/{name}")
    pd.DataFrame(res).T.sort_index().to_excel(f"{outdir}/{name}/forcing.xlsx")


def multi_shp(names, num_sample, tifs, shp_points_d, outdir):
    """计算多个 shapefile 的面均值"""
    for name in tqdm(names):
        one_shp(
            name,
            num_sample=num_sample,
            tifs=tifs,
            shp_points_d=shp_points_d,
            outdir=outdir,
        )


def main():
    folder_shp = "./folder_shp"
    folder_raster = "./folder_raster"
    outdir = "./output"

    shps = [x for x in absoluteFilePaths(folder_shp) if x.endswith(".shp")]
    tifs = absoluteFilePaths(folder_raster)
    shp_points_d = {}
    for shp in shps:
        name = shp.split("_")[-1].split(".")[0]
        points = list(np.round(shp_points(shp), 1))
        shp_points_d[name] = points
    names = list(shp_points_d.keys())

    proc = []

    num_threads = 8
    num_sample = 100000
    for i in range(num_threads):
        s, e = (len(names) // num_threads + 1) * i, (len(names) // num_threads + 1) * (
            i + 1
        )
        batch_names = names[s:e]
        p = Process(
            target=multi_shp, args=(batch_names, num_sample, tifs, shp_points_d, outdir)
        )
        proc.append(p)

    for p in proc:
        p.start()

    for p in proc:
        p.join()


if __name__ == "__main__":
    main()
