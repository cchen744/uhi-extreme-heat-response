# UHI Cell-Level Pipelin - Progress & Issues
## 1. GEE Spatial Filtering Performance (Partially resolved)
** Problem **
Filtering grid cells using:
```python
grid.filterBounds(region_proj)
```
was extremely slow and sometimes caused GEE to terminate the task
** Cause **
- grid contained thousands of polygon cells
- region_proj was a highly detailed city boundary polygon
- filterBounds() required many polygon-polygon intersection checks, which is computationally expensive in GEE
** Fix **
Simplify the city boundary before filtering
```python
region_proj = urban_region.transform(crs,1).simplify(500)
```
Also used centroid geometry in vector outputs:
```python
geometryType = "centroid"
```
This reduces geometry complexity and makes spatial filtering lighter
** Status **
Pipeline runs successfully with simplified geometry

## 2. incorrect Property Names from reduceToVectors() (Resolved)
** Problem **
Urban cell stats were missing from the final output table.
Expected fields:
- LST_urb_cells
- cell_n
- delta_uhi
But the output only contained:
```code
month
cell_id
LST_rur
rural_n
```
** Cause **
reduceToVectors() outputsreducer fields as:
```code
mean
count
```
instead of:
```code
{LST_band}_mean
{LST_band}_count
```
The code attempted to read incorrectly property names:
```python
ft.get(f"{LST_band}_mean")
ft.get(f"{LST_band}_count")
```
** Fix **
Updated property access to:
```python
lst_urb = ft.get("mean")
cell_n = ft.get("count")
```
** Status **
Urban cells stats now appear correctly in the final table

## 3. Interpretation of cell_n
It represents the # of valid MODIS pixels used to compute the cell_level LST.

## 4. Monthly Aggregation Strategy
Current approach:
```code
daily MODIS images
        ↓
QC filtering (clean_lst)
        ↓
monthly mean image
        ↓
cell-level aggregation
```
Instead of computing statistics for every daily image, the pipeline collapses daily LST images into a single monthly mean image
``` python
monthly_mean = ImageCollection(...)map(clean_lst).mean
```
This reduces computational cost
** Status **
Working as intended

## 5. Conceptual workflow (current pipeline)
1. Select city boundary
2. Build LCZ-based urban/rural masks
3. Compute monthly mean MODIS LST
4. Estimate rural reference temperature (city-level)
5. Compute urban LST for each 1km cell
6. Calculate:
```python
delta_uhi = LST_urb_cell - LST_rur
```
Note: this is still relatively slow - for the actual feature collection output, 1-month data costs 3 min's computation.

## Some Questions
1. Is the 3-12km rural ring an appropriate rural reference definition?
2. Is representing cells with centroid points instead of polygons acceptable?
3. Are there more efficient Earth Engine patterns for cell-level aggregations?

