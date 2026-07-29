"""
Reusable SUHI (Surface Urban Heat Island) Data Pipeline
Using Google Earth Engine (GEE) + PRISM (air temp) + Landsat 8/9 (LST)

Core Workflow:
1. Define urban region (city geometry) and rural reference ring (buffer annulus)
2. Label extreme-heat days using PRISM daily max air temperature — the
   90th percentile of summer (JJA) tmax, computed independently of any
   satellite LST product to avoid circularity.
3. For each Landsat 8/9 scene whose acquisition date falls in the extreme
   or baseline date list:
   - Mask cloud/shadow/cirrus pixels via QA_PIXEL
   - Convert ST_B10 to °C
   - Subtract that scene's own rural-ring mean temperature to get a
     per-scene UHI field (this makes scenes from different heat events
     comparable — otherwise averaging raw LST would reflect event
     intensity, not spatial pattern)
4. Average all per-scene UHI fields within each group (extreme vs. baseline)
   to get a 30m composite UHI map for each condition.
5. delta_UHI = mean(UHI | extreme) − mean(UHI | baseline)
"""
import ee
import geemap
import pandas as pd
import numpy as np
from datetime import datetime
from dateutil.relativedelta import relativedelta
from pathlib import Path

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
# 3. Build urban / rural masks (UA + LCZ)
# ------------------------------------------------------------------
def build_masks(city_geom, ring_outer_m, ring_inner_m, lcz_scale_m=100):

    urban_region = city_geom
    outer = city_geom.buffer(ring_outer_m)
    inner = city_geom.buffer(ring_inner_m)
    rural_region = outer.difference(inner)

    return urban_region, rural_region

# ------------------------------------------------------------------
# 5. Urban Area Selector
# ------------------------------------------------------------------

def select_ua(ua_fc, ua_contains=None, ua_name=None, ua_names=None):
    if ua_names is not None:
        return ua_fc.filter(ee.Filter.inList('NAME20', list(ua_names)))
    elif ua_name is not None:
        return ua_fc.filter(ee.Filter.eq('NAME20', ua_name))
    else:
        return ua_fc.filter(ee.Filter.stringContains('NAME20', ua_contains))


# ------------------------------------------------------------------
# 6. city-level wrapper to get PRISM + Landsat dataframe
# ------------------------------------------------------------------
def build_city_landsat_uhi(
    city_name,
    ua_fc,
    data_dir,
    start_date="2020-06-01",
    end_date="2026-07-28",
    ring_outer_m=12000,
    ring_inner_m=3000,
    extreme_percentile=90,
    cloud_max=40,
    export_scale_m=30,
    drive_folder="UHI_poster",
):
    """
    Runs the full PRISM-label + Landsat-composite pipeline for one city.

    Returns a summary dict (city, threshold, n_extreme_days, n_baseline_days,
    n_extreme_scenes, n_baseline_scenes, csv_path, export_task_names) and
    also writes:
      - {data_dir}/{city}_prism_labels.csv   — one row per PRISM day
      - two Export.image.toDrive tasks (started, not waited on) for the
        extreme and baseline UHI composites

    city_name can be a single NAME20 string or a list of NAME20 strings
    (for cities split across multiple TIGER Urban Area features, e.g. Phoenix).
    """
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)

    # 1) resolve city geometry
    if isinstance(city_name, (list, tuple)):
        city_fc = select_ua(ua_fc, ua_names=list(city_name))
        name_for_file = "_".join(city_name)
    else:
        city_fc = select_ua(ua_fc, ua_name=city_name)
        name_for_file = city_name
    city_geom = city_fc.geometry().dissolve()

    urban_region, rural_region = build_masks(city_geom, ring_outer_m, ring_inner_m)

    # 2) PRISM extreme-day labels
    prism_df, threshold = get_extreme_day_labels(
        city_geom, start_date, end_date, extreme_percentile=extreme_percentile
    )
    extreme_dates, baseline_dates = split_extreme_baseline_dates(prism_df)

    # save the label table — same role the old daily_suhi CSV played
    csv_path = data_dir / f"{name_for_file.replace(' ', '_').replace(',', '')}_prism_labels.csv"
    prism_df.to_csv(csv_path, index=False)
    print(f"{city_name}: saved PRISM labels to {csv_path}")

    # 3) Landsat UHI composites
    uhi_extreme, n_ext = landsat_uhi_composite(city_geom, rural_region, extreme_dates, cloud_max=cloud_max)
    uhi_base, n_base = landsat_uhi_composite(city_geom, rural_region, baseline_dates, cloud_max=cloud_max)

    n_ext_val = n_ext.getInfo()
    n_base_val = n_base.getInfo()
    print(f"{city_name}: Landsat scenes — extreme = {n_ext_val}, baseline = {n_base_val}")

    region = city_geom.buffer(3000)

    delta_uhi = uhi_extreme.subtract(uhi_base).rename("delta_UHI").toFloat()
    # start exporting GeoTIFF into GEE (visualization)
    task_ext = ee.batch.Export.image.toDrive(
        image=uhi_extreme.clip(region),
        description=f"{name_for_file}_UHI_extreme_{export_scale_m}m",
        folder=drive_folder,
        region=region,
        scale=export_scale_m,
        crs="EPSG:3857",
        maxPixels=1e9,
    )
    task_base = ee.batch.Export.image.toDrive(
        image=uhi_base.clip(region),
        description=f"{name_for_file}_UHI_baseline_{export_scale_m}m",
        folder=drive_folder,
        region=region,
        scale=export_scale_m,
        crs="EPSG:3857",
        maxPixels=1e9,
    )

    task_delta = ee.batch.Export.image.toDrive(
        image=delta_uhi.clip(region),
        description=f"{name_for_file}_UHI_delta_{export_scale_m}m",
        folder=drive_folder,
        region=region,
        scale=export_scale_m,
        crs="EPSG:3857",
        maxPixels=1e9,
    )
    task_ext.start()
    task_base.start()
    task_delta.start()

    print(f"{city_name}: export tasks started — check the Tasks tab in the GEE Code Editor "
          f"or run task.status() to monitor")

    return {
        "city": city_name,
        "threshold_C": threshold,
        "n_extreme_days": len(extreme_dates),
        "n_baseline_days": len(baseline_dates),
        "n_extreme_scenes": n_ext_val,
        "n_baseline_scenes": n_base_val,
        "csv_path": str(csv_path),
        "export_tasks": [task_ext.id, task_base.id, task_delta.id],
    }

# ------------------------------------------------------------------
# 7a. Extreme-day labels from PRISM only (no MODIS needed)
# ------------------------------------------------------------------
def get_extreme_day_labels(
    city_geom,
    start_date,
    end_date,
    extreme_percentile=90,
    summer_months=(6, 7, 8),
):
    """
    Returns a DataFrame: date | tmax | is_extreme, for every PRISM day in range.
 
    The threshold is the given percentile of tmax across summer days only,
    so 'extreme' means 'hot relative to this city's own summer', not
    relative to the whole year.
 
    PRISM is 2m air temperature interpolated from weather stations —
    completely independent of any satellite LST product, so the label
    is exogenous to the UHI signal being measured.
    """
    prism = ee.ImageCollection("OREGONSTATE/PRISM/ANd").select("tmax")
 
    def daily_mean(img):
        val = img.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=city_geom,
            scale=4000, # PRISM's initial resolution is 4km. Do 4000m to match that.
            maxPixels=1e9,
        ).get("tmax")
        return ee.Feature(None, {
            "date": img.date().format("YYYY-MM-dd"),
            "tmax": val,
        })
 
    fc = ee.FeatureCollection(
        prism.filterBounds(city_geom)
             .filterDate(start_date, end_date)
             .map(daily_mean)
    )
 
    df = geemap.ee_to_df(fc)
    df["date"] = pd.to_datetime(df["date"])
    df["tmax"] = pd.to_numeric(df["tmax"], errors="coerce") # If getting error, use 'NaN' as default filler
    df = df.dropna(subset=["tmax"]).sort_values("date")
 
    summer = df[df["date"].dt.month.isin(summer_months)]
    threshold = summer["tmax"].quantile(extreme_percentile / 100.0)
    df["is_extreme"] = ((df["tmax"] >= threshold) &
                        (df["date"].dt.month.isin(summer_months))).astype(int)
 
    n_extreme = int(df["is_extreme"].sum())
    n_summer = len(summer)
    print(f"PRISM p{extreme_percentile} threshold = {threshold:.2f} degC "
          f"(summer-only) | extreme days = {n_extreme} of {n_summer} summer days")
 
    return df, threshold

def split_extreme_baseline_dates(prism_df, summer_months=(6, 7, 8)):
    """
    Splits the PRISM label table into two explicit date lists for
    Landsat scene selection. Baseline = summer days that are NOT extreme.
    """
    summer = prism_df[prism_df["date"].dt.month.isin(summer_months)]
 
    extreme_dates = (summer.loc[summer["is_extreme"] == 1, "date"]
                           .dt.strftime("%Y-%m-%d").tolist())
    baseline_dates = (summer.loc[summer["is_extreme"] == 0, "date"]
                            .dt.strftime("%Y-%m-%d").tolist())
 
    print(f"extreme dates: {len(extreme_dates)} | baseline dates: {len(baseline_dates)}")
    return extreme_dates, baseline_dates

# ------------------------------------------------------------------
# 7b. Landsat UHI composite — per-scene normalization, then group mean
# ------------------------------------------------------------------
def landsat_uhi_composite(
    city_geom,
    rural_region,
    dates,
    cloud_max=40, # get the image if cloud coverage is smaller than 40%
    rural_scale_m=300,
):
    """
    Builds a mean UHI field (30m) from all Landsat 8/9 scenes whose
    acquisition date is in `dates`.
 
    Each scene is normalized against ITS OWN rural reference temperature
    before averaging. This matters: raw LST differs by several degrees
    between one heat event and another, so averaging raw LST across
    events would measure event intensity rather than spatial pattern.
    Subtracting each scene's own rural mean first makes scenes comparable.
 
    Returns (mean_uhi_image, scene_count).
    """
    landsat = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
               .merge(ee.ImageCollection("LANDSAT/LC09/C02/T1_L2"))) # Use two satellites ensure we are getting complete data
 
    def to_lst(img):
        qa = img.select("QA_PIXEL")
        clear = (qa.bitwiseAnd(1 << 3).eq(0)
                 .And(qa.bitwiseAnd(1 << 4).eq(0))
                 .And(qa.bitwiseAnd(1 << 2).eq(0)))
        lst = (img.select("ST_B10")
                  .multiply(0.00341802).add(149.0).subtract(273.15)
                  .rename("LST_C")
                  .updateMask(clear)) # convert raw digital number to Celsius
        return lst.copyProperties(img, ["system:time_start"])
 
    ic = (landsat
          .filterBounds(city_geom)
          .filter(ee.Filter.lt("CLOUD_COVER", cloud_max))
          .map(lambda i: i.set("date_str", i.date().format("YYYY-MM-dd")))
          .filter(ee.Filter.inList("date_str", ee.List(dates)))
          .map(to_lst))
 
    def tag_rural_ref(img):
        # coarser scale here is fine — this is a single scalar per scene,
        # and running it at 30m is needlessly slow
        ref = img.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=rural_region,
            scale=rural_scale_m,
            maxPixels=1e9, # Safety limit for reduceRegion() to make sure GEE doesn't crash
        ).get("LST_C")
        return img.set("rur_ref", ref)
 
    ic = ic.map(tag_rural_ref).filter(ee.Filter.notNull(["rur_ref"]))
 
    def to_uhi(img):
        return (img.subtract(ee.Number(img.get("rur_ref")))
                   .rename("UHI")
                   .toFloat()
                   .copyProperties(img, ["system:time_start"]))
 
    uhi_ic = ic.map(to_uhi)
    return uhi_ic.mean(), ic.size()