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

# ------------------------------------------------------------------
# 4. Daily aggregation (UA x day) - OPTIMIZED
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
# * 5. NEW FUNCTION: make_monthly_table_cells
# * Replaces make_grid_fc_2 + make_daily_table_cells entirely.
# *
# * Key design decisions:
# *   1. Monthly mean is computed in GEE (ImageCollection.mean()) before
# *      vectorizing, so reduceToVectors runs once per month, not once per day.
# *   2. Grid generation (pixelCoordinates → reduceToVectors) is inlined here
# *      and applied directly to the masked monthly mean image, eliminating
# *      the separate make_grid_fc_2 → reduceRegions two-step.
# *   3. delta_uhi = LST_urb_cell - LST_rur is computed server-side before
# *      export, so the output table is immediately analysis-ready.
# *   4. Output shape: (n_months × n_cells) instead of (n_days × n_cells),
# *      e.g. ~7k rows vs ~220k rows for Phoenix/1yr/500m.
# ------------------------------------------------------------------

def make_monthly_table_cells(
    start_date, end_date,
    urban_region, rural_region,
    urban_mask, rural_mask,
    lst_band, qc_band,
    lst_scale_m=1000,
    cell_scale_m=1000,       # * grid resolution in meters
    crs="EPSG:3857",
    err_m=100,
    tileScale=4
):
    """
    Returns a FeatureCollection with one feature per (month, cell):
        month        | YYYY-MM string
        cell_id      | integer pixel coordinate hash (x * 1e8 + y), stable across months
        LST_urb_cell | monthly mean urban LST for this cell (°C)
        cell_n       | number of valid MODIS pixels averaged into LST_urb_cell
        LST_rur      | monthly mean rural LST for the city (°C), one value per month
        rural_n      | number of valid MODIS pixels averaged into LST_rur
        delta_uhi    | LST_urb_cell - LST_rur (°C)
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
    region_proj = urban_region.transform(crs, 1)

    # * Combined reducer: mean + count in a single pass
    combined = ee.Reducer.mean().combine(ee.Reducer.count(), sharedInputs=True)

    monthly_fcs = []  # * collect one FC per month, flatten at the end

    for s, e in month_starts(start_date, end_date):
        month_str = s[:7]   # * "YYYY-MM"

        # * Step 1: collapse all valid daily images in this month to a single mean image
        monthly_mean = (
            ee.ImageCollection("MODIS/061/MYD11A1")
            .filterBounds(urban_region)
            .filterDate(s, e)
            .select([lst_band, qc_band])
            .map(lambda img: clean_lst(img, lst_band, qc_band))
            .mean()   # * temporal mean across days — reduces compute vs per-day loop
        )

        urb = monthly_mean.updateMask(urban_mask)
        rur = monthly_mean.updateMask(rural_mask)

        # * Step 2: rural reference — one scalar per month (cheap reduceRegion)
        rur_stats = rur.reduceRegion(
            reducer=combined,
            geometry=rural_region,
            scale=lst_scale_m,
            maxPixels=1e13
        )
        lst_rur = rur_stats.get(f"{lst_band}_mean")
        rural_n = rur_stats.get(f"{lst_band}_count")

        # * Step 3: urban cells — reduceToVectors on the monthly mean image.
        # *   Inlines the old make_grid_fc_2 + reduceRegions two-step into one call.
        # *   Running on a single monthly mean image is much faster than running
        # *   reduceToVectors on every day individually.
        urb_cells = (
            cell_id_img.addBands(urb)
               .reduceToVectors(
                   geometry=bounds,
                   scale=cell_scale_m,
                   geometryType="polygon",
                   crs=crs,
                   labelProperty="cell_id",
                   reducer=combined,
                   maxPixels=1e13,
                   tileScale=tileScale
               )
               .filterBounds(region_proj)
        )

        # * Step 4: attach month, rural reference, and delta_uhi to each cell feature
        def add_props(ft, _lst_rur=lst_rur, _rural_n=rural_n, _month=month_str):
            lst_urb  = ft.get(f"{lst_band}_mean")
            # * delta_uhi computed server-side so the exported CSV is analysis-ready
            delta    = ee.Number(lst_urb).subtract(ee.Number(_lst_rur))
            return ft.set({
                "month":        _month,
                "LST_urb_cell": lst_urb,
                "cell_n":       ft.get(f"{lst_band}_count"),
                "LST_rur":      _lst_rur,
                "rural_n":      _rural_n,
                "delta_uhi":    delta,
            }).select(["month", "cell_id", "LST_urb_cell", "cell_n",
                       "LST_rur", "rural_n", "delta_uhi"])

        monthly_fcs.append(urb_cells.map(add_props))

    # * Flatten list of monthly FCs into one FeatureCollection
    return ee.FeatureCollection(monthly_fcs).flatten()

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
    agg_func="median", # NOTE: only used for unit='city'; monthly cell path always uses mean
    ring_outer_m=12000,
    ring_inner_m=3000,
    lst_scale_m=1000,
    min_urban_pixels=50,
    min_rural_pixels=50,
    unit='city',
    cell_scale_m=1000,
    min_cell_pixels=1,
    cell_crs="EPSG:3857",
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
                # * MODIFIED: replaced make_daily_table_cells (day × cell) with
                # *           make_monthly_table_cells (month × cell).
                # *           The entire date range is handled in one call;
                # *           monthly iteration is done inside make_monthly_table_cells.
                fc = make_monthly_table_cells(
                    start_date, end_date,
                    urban_region, rural_region,
                    urban_mask, rural_mask,
                    lst_band, qc_band,
                    lst_scale_m=lst_scale_m,
                    cell_scale_m=cell_scale_m,
                    crs=cell_crs
                )
                fc = ee.FeatureCollection(fc)
                fc = fc.filter(ee.Filter.notNull(["month", "cell_id", "LST_urb_cell",
                                                    "cell_n", "LST_rur", "rural_n", "delta_uhi"]))
                fc = fc.filter(ee.Filter.gte("cell_n", min_cell_px))
                fc = fc.filter(ee.Filter.gte("rural_n", min_rural_px))
                fc = fc.select(["month", "cell_id", "LST_urb_cell", "cell_n",
                                  "LST_rur", "rural_n", "delta_uhi"])
                # * wrap in list to keep the rest of run_city's flattening logic intact
                fc_list_local.append(fc)
                month_ranges_local.append((start_date, end_date))

            else:
                raise ValueError("unit must be 'city' or 'cell'")
            
          return fc_list_local, month_ranges_local

    fc_list, month_ranges = build_fc_list(min_urban_pixels, min_rural_pixels, min_cell_pixels)  
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
        # * MODIFIED: dedup on (month, cell_id) instead of (date, cell_id);
        # *           delta_uhi already computed in GEE, no need to recalculate here
        if "month" not in df_all.columns:
            print("Missing 'month' column. Available:", list(df_all.columns))
            return pd.DataFrame()
        df_all = df_all.drop_duplicates(subset=["month", "cell_id"]).reset_index(drop=True)
        # delta_uhi is already in the DataFrame from GEE computation

    return df_all
