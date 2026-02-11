"""
Reusable SUHI data pipeline (GEE + MODIS LST)
Optimized Version: Uses Combined Reducers to reduce server load.
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
    lcz_img = ee.ImageCollection("RUB/RUBCLIM/LCZ/global_lcz_map/latest").first()
    lcz = lcz_img.select("LCZ_Filter")

    BUILT_MIN, BUILT_MAX = 1, 10
    WATER_CODE = 17

    is_built = lcz.gte(BUILT_MIN).And(lcz.lte(BUILT_MAX))
    is_water = lcz.eq(WATER_CODE)
    is_natural = is_built.Not().And(is_water.Not())

    urban_region = city_geom
    outer = city_geom.buffer(ring_outer_m)
    inner = city_geom.buffer(ring_inner_m)
    rural_region = outer.difference(inner)

    urban_mask = is_water.Not().clip(urban_region)
    rural_mask = is_natural.clip(rural_region)

    return urban_region, rural_region, urban_mask, rural_mask

# ================ Here we use ReduceToVector to simplify the grid cell geenration ============
def make_grid_fc_2(region_geom, cell_size_m=1000, crs="EPSG:3857", err_m=100, tileScale=4):
    """
    Build 1km x 1km fishnet polygons in a metric CRS using reduceToVectors.
    Much more stable than List.sequence + map + flatten.
    """
    # 1) Bounds in target CRS (meters)
    bounds = region_geom.bounds(ee.ErrorMargin(err_m)).transform(crs, 1)
    region_proj = region_geom.transform(crs, 1)

    # 2) Create a constant image in the desired projection/scale
    proj = ee.Projection(crs).atScale(cell_size_m)
    pc = ee.Image.pixelCoordinates(proj)
    x = pc.select("x").toInt()
    y = pc.select("y").toInt()

    # Build a per-pixel unique-ish id so adjacent pixels won't merge.
    # (Use a big multiplier; safe for typical WebMercator meter ranges.)
    pid = x.multiply(100000000).add(y).rename("pid")

    #vadd a dummy value band so Reducer.first has something to reduce
    img2 = pid.addBands(ee.Image.constant(1).rename("v"))

    # 3) Vectorize pixels inside bounds -> grid polygons
    grid = img2.reduceToVectors(
        geometry=bounds,
        scale=cell_size_m,
        geometryType="polygon",
        crs=crs,
        labelProperty="pid", 
        reducer=ee.Reducer.first(),
        maxPixels=1e13,
        tileScale=tileScale
    )

    
    grid = grid.filterBounds(region_proj)


    # # 4) Clip each cell to the city and create a stable-ish cell_id from centroid coords
    # def post(ft):
    #     ft = ee.Feature(ft)
    #     geom = ft.geometry()
    #     inter = geom.intersects(region_proj, ee.ErrorMargin(err_m))
    #     geom_clip = ee.Geometry(
    #         ee.Algorithms.If(
    #             inter,
    #             geom.intersection(region_proj, ee.ErrorMargin(err_m)),
    #             geom
    #         )
    #     )

    #     # centroid in CRS -> use rounded meters as id components
    #     c = geom.centroid(ee.ErrorMargin(err_m)).transform(crs, 1).coordinates()
    #     x = ee.Number(ee.List(c).get(0)).round()
    #     y = ee.Number(ee.List(c).get(1)).round()
    #     cell_id = x.format('%.0f').cat('_').cat(y.format('%.0f'))

    #     return ee.Feature(geom_clip, {"cell_id": cell_id}).set("keep", inter)
    
    # give every cell a stable cell_id using its centroid in "EPSG:3857"
    def add_id(ft):
      ft = ee.Feature(ft)
      c = ft.geometry().centroid(ee.ErrorMargin(err_m)).transform(crs, 1).coordinates()
      x = ee.Number(ee.List(c).get(0)).round()
      y = ee.Number(ee.List(c).get(1)).round()
      cell_id = x.format('%.0f').cat('_').cat(y.format('%.0f'))
      return ft.set({"cell_id": cell_id})

    return grid.map(add_id).select(["cell_id"])


# ------------------------------------------------------------------
# 5. Daily aggregation (UA x day) - OPTIMIZED
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

    # OPTIMIZATION#1: Combine Mean/Median with Count into a single Reducer
    # This reduces the number of reduceRegion calls by half.
    base_reducer = ee.Reducer.mean() if agg_func == "mean" else ee.Reducer.median()
    combined_reducer = base_reducer.combine(
        reducer2=ee.Reducer.count(),
        sharedInputs=True
    )

    def agg(img):
        date = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd")

        urb = img.updateMask(urban_mask)
        rur = img.updateMask(rural_mask)

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

    return ic.map(agg)

# ------------------------------------------------------------------
# 6. UPDATE: new function for grid-cell-level daily output
# ------------------------------------------------------------------
def make_daily_table_cells(
    start, end,
    urban_region, rural_region,
    urban_mask, rural_mask,
    lst_band, qc_band,
    agg_func="mean",
    lst_scale_m=1000,
    cell_scale_m=1000,
    crs="EPSG:3857"
):
    """
    Output FeatureCollection with one feature per (day, cell).
    Each feature contains:
      - date
      - cell_id
      - LST_urb_cell, cell_n
      - LST_rur (city-level), rural_n
    """

    # 1) Build 1km fishnet grid within urban_region (projected to CRS meters)
    grid_fc = make_grid_fc_2(urban_region, cell_size_m=cell_scale_m, crs=crs)

    # 2) Load & clean MODIS LST images for the given time window
    ic = (
        ee.ImageCollection("MODIS/061/MYD11A1")
        .filterBounds(urban_region)
        .filterDate(start, end)
        .select([lst_band, qc_band])
        .map(lambda img: clean_lst(img, lst_band, qc_band))
    )

    # 3) Combined reducer: mean/median + count in one pass (same as city-level)
    base_reducer = ee.Reducer.mean() if agg_func == "mean" else ee.Reducer.median()
    combined_reducer = base_reducer.combine(
        reducer2=ee.Reducer.count(),
        sharedInputs=True
    )

    key_val = agg_func      # "mean" or "median"
    key_cnt = "count"

    # new：rural reduceRegion
    rur_key_val_band = f"{lst_band}_{agg_func}"   # e.g. LST_Night_1km_mean
    rur_key_cnt_band = f"{lst_band}_count"        # e.g. LST_Night_1km_count

    def agg(img):
        # 4) Attach a date string to every output row (per day)
        date = ee.Date(img.get("system:time_start")).format("YYYY-MM-dd")

        # 5) Apply masks: urban excludes water; rural uses natural-only ring
        urb = img.updateMask(urban_mask)
        rur = img.updateMask(rural_mask)

        # 6) Rural reference is city-level (one value per day)
        rur_stats = rur.reduceRegion(
            reducer=combined_reducer,
            geometry=rural_region,
            scale=lst_scale_m,
            maxPixels=1e13
        )
        
        lst_rur = ee.Algorithms.If(
          rur_stats.contains(rur_key_val_band),
          rur_stats.get(rur_key_val_band),
          ee.Algorithms.If(rur_stats.contains(key_val), rur_stats.get(key_val), None)
        )

        rural_n = ee.Algorithms.If(
            rur_stats.contains(rur_key_cnt_band),
            rur_stats.get(rur_key_cnt_band),
            ee.Algorithms.If(rur_stats.contains(key_cnt), rur_stats.get(key_cnt), 0)
        )

        # 7) Urban is reduced per grid cell (many values per day)
        #    reduceRegions returns a FeatureCollection with stats per feature geometry.
        urb_cells = urb.reduceRegions(
            collection=grid_fc,
            reducer=combined_reducer,
            scale=lst_scale_m,
            tileScale=4
        )

        # 8) For each cell, write city-level rural reference into its properties
        def add_rural_props(ft):
            ft = ee.Feature(ft)
            return ft.set({
                "date": date,
                "LST_rur": lst_rur,
                "rural_n": rural_n,
                "LST_urb_cell": ft.get(key_val),
                "cell_n": ft.get(key_cnt),
            })

        urb_cells = urb_cells.map(add_rural_props)

        # 9) Keep only needed columns
        return urb_cells.select(["date", "cell_id", "LST_urb_cell", "cell_n", "LST_rur", "rural_n"])

    # 10) Map over days and flatten day-wise FC into one FC
    return ee.FeatureCollection(ic.map(agg)).flatten()


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

# ------------------------------------------------------------------
# 7. UPDATE: Main entry: have options for running one city or running one grid cell
# ------------------------------------------------------------------
def run_city(
    *,
    ua_fc,
    start_date,
    end_date,
    ua_name=None,
    ua_contains=None,
    ua_names=None,
    lst_band="LST_Night_1km",
    qc_band="QC_Night",
    agg_func="median",
    ring_outer_m=12000,
    ring_inner_m=3000,
    lst_scale_m=1000,
    min_urban_pixels=50,
    min_rural_pixels=50,
    unit='city',
    cell_scale_m=1000,
    min_cell_pixels=1,
    cell_crs="EPSG:3857",
    extreme_percentile=90,
    out_csv=None,
    export_to_drive=False,
    export_desc=None,
    export_folder="UHI_exports",
    debug=False,
    relax_filters_on_empty=True
):
    city_fc = select_ua(ua_fc, ua_name=ua_name, ua_contains=ua_contains, ua_names=ua_names)
    if city_fc.size().getInfo() == 0:
        raise ValueError("No UA matched your query.")
    city_geom = city_fc.geometry()

    urban_region, rural_region, urban_mask, rural_mask = build_masks(
        city_geom, ring_outer_m, ring_inner_m
    )

    def build_fc_list(min_urban_px, min_rural_px, min_cell_px):
          fc_list_local = []
          month_ranges_local = []
          for s, e in month_starts(start_date, end_date):
            month_ranges_local.append((s, e))
            if unit == 'city':
              fc = make_daily_table(
                  s, e,
                  urban_region, rural_region,
                  urban_mask, rural_mask,
                  lst_band, qc_band,
                  agg_func, lst_scale_m
              )

              fc = ee.FeatureCollection(fc)
              # Clean up nulls
              fc = fc.filter(ee.Filter.notNull(["LST_urb", "LST_rur", "urban_n", "rural_n", "date"]))
              fc = fc.filter(ee.Filter.gte("urban_n", min_urban_px))
              fc = fc.filter(ee.Filter.gte("rural_n", min_rural_px))
              fc = fc.select(["date", "LST_urb", "LST_rur", "urban_n", "rural_n"])
            
            elif unit == "cell":
              fc = make_daily_table_cells(
              s, e,
              urban_region, rural_region,
              urban_mask, rural_mask,
              lst_band, qc_band,
              agg_func,
              lst_scale_m=lst_scale_m,
              cell_scale_m=cell_scale_m,
              crs=cell_crs
              )
              fc = ee.FeatureCollection(fc)
              fc = fc.filter(ee.Filter.notNull(["date", "cell_id", "LST_urb_cell", "cell_n", "LST_rur", "rural_n"]))
              fc = fc.filter(ee.Filter.gte("cell_n", min_cell_px))
              fc = fc.filter(ee.Filter.gte("rural_n", min_rural_px))
              fc = fc.select(["date", "cell_id", "LST_urb_cell", "cell_n", "LST_rur", "rural_n"])


            else:
              raise ValueError("unit must be 'city' or 'cell'")
              
            fc_list_local.append(fc)

            return fc_list_local, month_ranges_local


    fc_list, month_ranges = build_fc_list(min_urban_pixels, min_rural_pixels, min_cell_pixels)  
    # checkpoint 1
    print("fc_list length:", len(fc_list))
    print("first fc size:", fc_list[0].size().getInfo())     

    if not fc_list:
      return pd.DataFrame()
  
    fc_all = ee.FeatureCollection(fc_list).flatten()
    fc_all_size = fc_all.size().getInfo()
    if fc_all_size == 0:
      if relax_filters_on_empty:
        relaxed_urban_px = 1
        relaxed_rural_px = 1
        relaxed_cell_px = 1
        print(
          "No features returned after filtering. Retrying with relaxed pixel thresholds: "
          f"urban>={relaxed_urban_px}, rural>={relaxed_rural_px}, cell>={relaxed_cell_px}."
        )
        fc_list, month_ranges = build_fc_list(relaxed_urban_px, relaxed_rural_px, relaxed_cell_px)
        fc_all = ee.FeatureCollection(fc_list).flatten()
        fc_all_size = fc_all.size().getInfo()
      if fc_all_size == 0:
        print("No features returned from Earth Engine after filtering. Returning empty DataFrame.")
        if debug:
          print("Monthly feature counts after filtering:")
          for (s, e), fc in zip(month_ranges, fc_list):
            monthly_size = fc.size().getInfo()
            print(f"  {s} to {e}: {monthly_size}")
          print(
            "Consider lowering min_urban_pixels/min_rural_pixels "
            "(or min_cell_pixels for unit='cell') or expanding date/AOI filters."
          )
        return pd.DataFrame()

    fc_all_size = fc_all.size().getInfo()
    if fc_all_size == 0:
      print("No features returned from Earth Engine after filtering. Returning empty DataFrame.")
      return pd.DataFrame()

    if export_to_drive:
      desc = export_desc or f"UHI_{ua_name or ua_contains or 'city'}"
      task = ee.batch.Export.table.toDrive(
        collection=fc_all,
        description=desc,
        fileFormat="CSV",
        fileNamePrefix=desc,
        folder=export_folder
      )
      task.start()
      return None

    # Optimized: Try to fetch data. 
    # The previous 500 error happens here. The upstream optimization (combined reducer) should help.
    try:
        df_all = geemap.ee_to_df(fc_all)
        print("df_all columns:", list(df_all.columns))
        print("df_all shape:", df_all.shape)
        print(df_all.empty)
    except Exception as e:
        print(f"Error fetching data: {e}")
        return pd.DataFrame()
    if "date" not in df_all.columns:
      print("Missing 'date' column in DataFrame. Available columns:", list(df_all.columns))
      return pd.DataFrame()

    df_all["date"] = pd.to_datetime(df_all["date"])


    if unit == "city":
      df_all = df_all.drop_duplicates(subset=["date"]).reset_index(drop=True)
      df_all["SUHI"] = df_all["LST_urb"] - df_all["LST_rur"]


    elif unit == "cell":
      df_all = df_all.drop_duplicates(subset=["date", "cell_id"]).reset_index(drop=True)
      df_all["SUHI"] = df_all["LST_urb_cell"] - df_all["LST_rur"]

    return df_all
