"""
Reusable SUHI data pipeline (GEE + MODIS LST)

Optimized Version: GEE Cell Monthly Aggregation --> Notebook Cell Monthly delta_SUHI

Core Workflow:
Each month → Apply ImageCollection.mean() to LST images for all days in that month → Obtain a monthly mean image
→ Perform reduceToVectors on the monthly mean image (aggregate by grid cells) → Obtain monthly mean urban LST for each cell
→ Simultaneously apply reduceRegion to rural regions → Obtain monthly mean rural LST
→ delta_uhi = urban_cell_mean - rural_mean
"""

import ee
import geemap
import pandas as pd
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta

# ------------------------------------------------------------------
# 0. Earth Engine init
# ------------------------------------------------------------------
def init_ee():
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
    
    # bits 0-1: mandatory QA
    mandatory = qc.bitwiseAnd(3)
    produced = mandatory.lte(1)
    valid = lst.neq(0)
    good = produced.And(valid)
    
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
    # lcz_img = ee.ImageCollection("RUB/RUBCLIM/LCZ/global_lcz_map/latest").first()
    # lcz = lcz_img.select("LCZ_Filter")

    # BUILT_MIN, BUILT_MAX = 1, 10
    # WATER_CODE = 17

    # is_built = lcz.gte(BUILT_MIN).And(lcz.lte(BUILT_MAX))
    # is_water = lcz.eq(WATER_CODE)
    # is_natural = is_built.Not().And(is_water.Not())

    urban_region = city_geom
    outer = city_geom.buffer(ring_outer_m)
    inner = city_geom.buffer(ring_inner_m)
    rural_region = outer.difference(inner)

    return urban_region, rural_region

# ------------------------------------------------------------------
# 4. Daily aggregation (UA x day) - OPTIMIZED
# ------------------------------------------------------------------
def make_daily_table(
    start, end,
    urban_region, rural_region,
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

    # OPTIMIZATION#1: Combine Mean/Median with Count into a single Reducer
    # This reduces the number of reduceRegion calls by half.
    base_reducer = ee.Reducer.mean() if agg_func == "mean" else ee.Reducer.median()
    combined_reducer = base_reducer.combine(
        reducer2=ee.Reducer.count(),
        sharedInputs=True
    )

    def agg(img):
        date = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd")

        urb = img.clip(urban_region)
        rur = img.clip(rural_region)

        # 1. Urban Reduction (Both Stats at once)
        urb_stats = urb.reduceRegion(
            reducer=combined_reducer,
            geometry=urban_region,
            scale=lst_scale_m,
            maxPixels=1e13
        )

        # 2. Rural Reduction (Both Stats at once)
        rur_stats = rur.reduceRegion(
            reducer=combined_reducer,
            geometry=rural_region,
            scale=lst_scale_m,
            maxPixels=1e13
        )
        
        # Construct Output Keys based on default reducer outputs
        key_val = f"{lst_band}_mean" if agg_func == "mean" else f"{lst_band}_median"
        key_cnt = f"{lst_band}_count"

        return ee.Feature(None, {
            "date": date,
            "LST_urb": urb_stats.get(key_val),
            "LST_rur": rur_stats.get(key_val),
            "urban_n": urb_stats.get(key_cnt),
            "rural_n": rur_stats.get(key_cnt),
        })

    final_fc = ee.FeatureCollection(ic.map(agg))
    return final_fc

def make_daily_table_cells(
    start_date, end_date,
    urban_region, rural_region,
    lst_band, qc_band,
    lst_scale_m=1000,
    cell_scale_m=1000,       # * grid resolution in meters
    crs="EPSG:3857",
    err_m=100,
    tileScale=4
):
    """
    Returns a FeatureCollection with one image per (date, cell):
        date        | YYYY-MM-dd string
        cell_id      | integer pixel coordinate hash (x * 1e8 + y), stable across days
        LST_urb_cell | daily mean urban LST for this cell (°C)
        urb_cell_n       | number of valid MODIS pixels averaged into LST_urb_cell
        LST_rur      | daily mean rural LST for the city (°C), one value one day
        rural_cell_n      | number of valid MODIS pixels averaged into LST_rur
        uhi    | LST_urb_cell - LST_rur (°C)
    """
    # * Precompute grid geometry constants (shared across all months)
    proj       = ee.Projection(crs).atScale(cell_scale_m)
    pc         = ee.Image.pixelCoordinates(proj)
    # * cell_id image: integer hash of pixel (x, y) coordinates — unique and stable
    cell_id_img = (
        pc.select("x").toInt()
          .multiply(100000000)
          .add(pc.select("y").toInt())
          .rename("cell_id")
    )
    bounds      = urban_region.bounds(ee.ErrorMargin(err_m)).transform(crs, 1)
    region_proj = urban_region.transform(crs, 1).simplify(cell_scale_m / 2)

    # * Combined reducer: mean + count in a single pass
    combined = ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True)
    
    dummy = ee.Image.constant(1).rename("dummy")
    grid = (
          cell_id_img.addBands(dummy)
          .reduceToVectors(
                  geometry=bounds,
                  scale=cell_scale_m,
                  geometryType="polygon", # Determining whether a point lies inside a polygon is more efficient than determining whether it intersects the polygon.
                  crs=crs,
                  labelProperty="cell_id",
                  reducer=ee.Reducer.first(),
                  maxPixels=1e13,
                  tileScale=tileScale
              )
              .filterBounds(region_proj)
      )

    # * Step 1: collect daily image
    daily_ics = (
        ee.ImageCollection("MODIS/061/MYD11A1")
        .filterBounds(urban_region)
        .filterDate(start_date, end_date)
        .select([lst_band, qc_band])
        .map(lambda img: clean_lst(img, lst_band, qc_band))
        # * temporal mean across days — reduces compute vs per-day loop
    )

    # print("Rural stat:", rur_stats.getInfo())
    
    daily_fc_list = [] # container of final output
    img_list = daily_ics.toList(daily_ics.size()) # a list of daily images

    def single_img_processer(img):
      '''
      for each image, generate a grided imge with uhi in each grid cell
      '''
      date_str = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd")
      urb = img.clip(urban_region)
      rur = img.clip(rural_region)

      # * Step 1: rural reference — one scalar per day
      rur_stats = rur.reduceRegion(
        reducer=combined,
        geometry=rural_region,
        scale=lst_scale_m,
        maxPixels=1e13
    )
      lst_rur = rur_stats.get(f"{lst_band}_mean")
      rural_cell_n = rur_stats.get(f"{lst_band}_count")

      # * Step 2: urban cells —  on the daily image.
      
    
      urb_cells = urb.reduceRegions(
          collection=grid,
          reducer=combined,
          scale=lst_scale_m,
          crs=img.projection(),
          tileScale=tileScale
      )

      # print("None empty urb cell:", urb_cells
     # .filter(ee.Filter.gt("count", 0))
     # .first()
     # .getInfo())

      # Step 3: attach daily, rural reference, and uhi to each cell feature
      def add_props(ft, _lst_rur=lst_rur, _rural_cell_n=rural_cell_n, _date=date_str):
          lst_urb = ft.get("mean")
          urb_cell_n = ft.get("count") # cell_n = The number of valid MODIS pixels used to compute the temperature for this cell
          # uhi computed server-side so the exported CSV is analysis-ready
          uhi = ee.Algorithms.If(
              ee.Algorithms.IsEqual(lst_urb, None),
              None,
              ee.Algorithms.If(
                  ee.Algorithms.IsEqual(_lst_rur, None),
                  None,
                  ee.Number(lst_urb).subtract(ee.Number(_lst_rur))
              )
          )
          return ft.set({
              "date":        _date,
              "LST_urb_cell": lst_urb,
              "urb_cell_n":       urb_cell_n,
              "LST_rur":      _lst_rur,
              "rural_cell_n":      _rural_cell_n,
              "uhi":    uhi,
          }).select(["date", "cell_id", "LST_urb_cell", "urb_cell_n",
                      "LST_rur", "rural_cell_n", "uhi"])
      
      mapped_fc = urb_cells.map(add_props)

      # print("Non-empty urb cells with properties:", mapped_fc
      # .filter(ee.Filter.gt("urb_cell_n", 0))
      # .first()
      # .getInfo())

      return mapped_fc

    n = daily_ics.size().getInfo()

    for i in range(n):
      img = ee.Image(img_list.get(i))
      daily_fc_list.append(single_img_processer(img))
    # Fixing empty dataset:filterDate is silently dropping everything because features use a string 'date' property instead of system:time_start.
    def add_time_start(feature):
        date = ee.Date(feature.get('date'))
        return feature.set('system:time_start', date.millis())
    final_fc = ee.FeatureCollection(daily_fc_list).flatten()
    
    return final_fc.map(add_time_start)

# ------------------------------------------------------------------
# 6. Urban Area Selector
# ------------------------------------------------------------------

def select_ua(ua_fc, *, ua_name=None, ua_contains=None, ua_names=None):
    args = [ua_name is not None, ua_contains is not None, ua_names is not None]
    if sum(args) != 1:
        raise ValueError("Provide exactly one of ua_name, ua_contains, ua_names")

    if ua_name is not None:
        fc = ua_fc.filter(ee.Filter.eq("NAME20", ua_name))
    elif ua_contains is not None:
        fc = ua_fc.filter(ee.Filter.stringContains("NAME20", ua_contains))
    else:
        filt = None
        for name in ua_names:
            f = ee.Filter.eq("NAME20", name)
            filt = f if filt is None else ee.Filter.Or(filt, f)
        fc = ua_fc.filter(filt)
    return fc
