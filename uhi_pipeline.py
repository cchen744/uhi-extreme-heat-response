"""
Reusable SUHI data pipeline (GEE + MODIS LST)
Extracted from 01_data_exploration.ipynb
"""

import ee
import geemap
import pandas as pd
import numpy as np
import os
from datetime import datetime
from dateutil.relativedelta import relativedelta


# ------------------------------------------------------------------
# 0. Earth Engine init
# ------------------------------------------------------------------
def init_ee():
    """Initialize Earth Engine (call once per runtime)."""
    try:
        ee.Initialize()
    except Exception:
        ee.Authenticate()
        ee.Initialize()


# ------------------------------------------------------------------
# 1. Helper: monthly ranges
# ------------------------------------------------------------------
def month_starts(start_date, end_date):
    s = datetime.strptime(start_date, "%Y-%m-%d").replace(day=1)
    e = datetime.strptime(end_date, "%Y-%m-%d")
    cur = s
    while cur < e:
        nxt = cur + relativedelta(months=1)
        yield cur.strftime("%Y-%m-%d"), nxt.strftime("%Y-%m-%d")
        cur = nxt


# ------------------------------------------------------------------
# 2. MODIS LST cleaning (QA)
# ------------------------------------------------------------------
def clean_lst(img, lst_band, qc_band):
    lst = img.select(lst_band)
    qc = img.select(qc_band).toUint16()

    # bits 0–1: mandatory QA
    mandatory = qc.bitwiseAnd(3)
    produced = mandatory.lte(1)   # keep 0 or 1
    valid = lst.neq(0)            # drop fill

    good = produced.And(valid)

    # scale: 0.02 K → °C
    lst_c = lst.multiply(0.02).subtract(273.15)

    return (
        lst_c.updateMask(good)
             .rename(lst_band)
             .copyProperties(img, ["system:time_start"])
    )


# ------------------------------------------------------------------
# 3. Build urban / rural masks (UA + LCZ)
# ------------------------------------------------------------------
def build_masks(city_geom, ring_outer_m, ring_inner_m, lcz_scale_m=100):
    """
    Urban:
      - inside UA
      - exclude water
    Rural:
      - ring outside UA
      - natural LCZ only (exclude built + water)
    """
    # LCZ global map
    lcz_img = ee.ImageCollection(
        "RUB/RUBCLIM/LCZ/global_lcz_map/latest"
    ).first()
    lcz = lcz_img.select("LCZ_Filter")

    BUILT_MIN, BUILT_MAX = 1, 10
    WATER_CODE = 17

    is_built = lcz.gte(BUILT_MIN).And(lcz.lte(BUILT_MAX))
    is_water = lcz.eq(WATER_CODE)
    is_natural = is_built.Not().And(is_water.Not())

    # regions
    urban_region = city_geom
    outer = city_geom.buffer(ring_outer_m)
    inner = city_geom.buffer(ring_inner_m)
    rural_region = outer.difference(inner)

    # masks
    urban_mask = is_water.Not().clip(urban_region)
    rural_mask = is_natural.clip(rural_region)

    return urban_region, rural_region, urban_mask, rural_mask


# ------------------------------------------------------------------
# 4. Daily aggregation (UA × day)
# ------------------------------------------------------------------
def make_daily_table(
    start, end,
    urban_region, rural_region,
    urban_mask, rural_mask,
    lst_band, qc_band,
    agg_func="mean",
    lst_scale_m=1000
):
    ic = (
        ee.ImageCollection("MODIS/061/MYD11A1")
        .filterBounds(urban_region)
        .filterDate(start, end)
        .select([lst_band, qc_band])
        .map(lambda img: clean_lst(img, lst_band, qc_band))
    )

    reducer = ee.Reducer.mean() if agg_func == "mean" else ee.Reducer.median()

    def agg(img):
        date = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd")

        urb = img.updateMask(urban_mask)
        rur = img.updateMask(rural_mask)

        urb_stat = urb.reduceRegion(
            reducer=reducer,
            geometry=urban_region,
            scale=lst_scale_m,
            maxPixels=1e13
        )

        rur_stat = rur.reduceRegion(
            reducer=reducer,
            geometry=rural_region,
            scale=lst_scale_m,
            maxPixels=1e13
        )

        urb_n = urb.reduceRegion(
            reducer=ee.Reducer.count(),
            geometry=urban_region,
            scale=lst_scale_m,
            maxPixels=1e13
        )

        rur_n = rur.reduceRegion(
            reducer=ee.Reducer.count(),
            geometry=rural_region,
            scale=lst_scale_m,
            maxPixels=1e13
        )

        return ee.Feature(None, {
            "date": date,
            "LST_urb": urb_stat.get(lst_band),
            "LST_rur": rur_stat.get(lst_band),
            "urban_n": urb_n.get(lst_band),
            "rural_n": rur_n.get(lst_band),
        })

    return ic.map(agg)


# ------------------------------------------------------------------
# 5. Main entry: run one city
# ------------------------------------------------------------------
def run_city(
    city_name,
    ua_fc,
    start_date,
    end_date,
    lst_band="LST_Night_1km",
    qc_band="QC_Night",
    agg_func="mean",
    ring_outer_m=30000,
    ring_inner_m=5000,
    lst_scale_m=1000,
    min_urban_pixels=30,
    min_rural_pixels=30,
    extreme_percentile=90,
    out_csv=None
):
    # get city geometry
    city = ua_fc.filter(ee.Filter.stringContains("NAME20", city_name))
    city_geom = city.geometry()

    # masks
    urban_region, rural_region, urban_mask, rural_mask = build_masks(
        city_geom,
        ring_outer_m,
        ring_inner_m
    )

    dfs = []

    for s, e in month_starts(start_date, end_date):
        fc = make_daily_table(
            s, e,
            urban_region, rural_region,
            urban_mask, rural_mask,
            lst_band, qc_band,
            agg_func, lst_scale_m
        )

        fc = ee.FeatureCollection(fc)
        df = geemap.ee_to_df(fc)

        if df.empty:
            continue

        df = df.dropna()
        df = df[
            (df["urban_n"] >= min_urban_pixels) &
            (df["rural_n"] >= min_rural_pixels)
        ]

        if not df.empty:
            dfs.append(df)

    if not dfs:
        return pd.DataFrame()

    df_all = pd.concat(dfs, ignore_index=True)
    df_all["date"] = pd.to_datetime(df_all["date"])
    df_all = (
        df_all
        .sort_values("date")
        .drop_duplicates(subset=["date"])
        .reset_index(drop=True)
    )

    df_all["SUHI"] = df_all["LST_urb"] - df_all["LST_rur"]

    # Extreme heat (percentile-based, no consecutive requirement)
    thr = np.percentile(df_all["LST_urb"], extreme_percentile)
    df_all["is_extreme"] = (df_all["LST_urb"] >= thr).astype(int)

    if out_csv:
        df_all.to_csv(out_csv, index=False)

    return df_all
