# Variable Definitions

This project investigates whether Urban Heat Island (UHI) intensity is amplified during extreme heat events, and whether such amplification is associated with built environment characteristics.

All variables used in the pipeline are defined below.

## 1. Spatial Units
**UAC (Urbanized Area Cluster)**

A polygon representing the urban footprint of Phoenix.
This defines the urban study area.

**Rural Reference Ring**

An area outside the urban footprint used to estimate background temperature.

It is defined as the region between an outer buffer and an inner buffer around the urban footprint.

Example configuration:

```code
Inner buffer: 3 km from the UAC boundary

Outer buffer: 12 km from the UAC boundary
```

The rural reference ring is the area between these two buffers.

**Urban Grid Cells**

The urban area (UAC) is divided into equal-sized grid cells.

Typical configuration:'

```code
Grid resolution ≈ 1 km, consistent with the spatial resolution of MODIS LST.
```

Each grid cell represents the spatial unit of analysis.

## 2. City-Level Temperature (Event Definition)
**T_city_day**

Daily mean air temperature over the urban system.

This is calculated as the spatial average of air temperature (for example from PRISM) across the UAC area.

Purpose:

Used only to define extreme heat events.

Not used directly in UHI calculation.

## 3. Extreme Heat Definition
**extreme_threshold**

A temperature threshold used to identify extreme heat days.

Example:

```code
The 90th percentile of the city-level temperature time series.
```

**Extreme Days**

Days when the city-level temperature is greater than or equal to the extreme threshold.

These represent periods of extreme heat conditions.

**Baseline Days**

Days that do not meet the extreme heat criterion.

These represent normal background conditions used for comparison.

## 4. Surface Temperature Variables

Surface temperature is derived from MODIS Land Surface Temperature (LST).

**T_cell_day(i)**

Daily land surface temperature for urban grid cell i.

Each cell receives the LST value of the MODIS pixel intersecting that grid cell.

**T_rural_day**

Daily mean surface temperature over the rural reference ring.

This is calculated as the spatial average of LST across the rural ring area.

## 5. Aggregated Temperature Variables

Temperatures are averaged separately for extreme days and baseline days.

**T_cell_extreme(i)**

Mean surface temperature of cell i during extreme heat days.

**T_cell_baseline(i)**

Mean surface temperature of cell i during baseline (non-extreme) days.

**T_rural_extreme**

Mean rural surface temperature during extreme heat days.

**T_rural_baseline**

Mean rural surface temperature during baseline days.

## 6. Urban Heat Island (UHI)
**Cell-Level UHI Intensity**

UHI intensity measures the temperature difference between an urban location and the rural background.

```code
For extreme conditions:
UHI_cell_extreme(i) = T_cell_extreme(i) − T_rural_extreme
For baseline conditions:
UHI_cell_baseline(i) = T_cell_baseline(i) − T_rural_baseline
```

## 7. UHI Intensification

The main outcome variable measures whether UHI becomes stronger during extreme heat events.

```code
ΔUHI_cell(i) = UHI_cell_extreme(i) − UHI_cell_baseline(i)
```

Interpretation:

ΔUHI_cell > 0 → UHI intensifies during extreme heat

ΔUHI_cell = 0 → no change in UHI intensity

ΔUHI_cell < 0 → UHI weakens during extreme heat

## 8. Built Environment Variables

Built environment characteristics are extracted for each urban grid cell.

Examples include:

```code
LCZ_type(i)
Local Climate Zone classification describing urban form and land cover.

impervious_fraction(i)
Percentage of impervious surface within the grid cell.

NDVI_mean(i)
Average vegetation index within the grid cell.

built_fraction(i)
Fraction of built-up land cover.

elevation(i)
Mean elevation within the grid cell.
```

These variables are used to explain spatial variation in ΔUHI_cell.

## 9. Final Analysis Dataset

Each row in the final dataset represents one urban grid cell.

Example structure:

cell_id | T_cell_extreme | T_cell_baseline | UHI_extreme	| UHI_baseline	| ΔUHI_cell	| LCZ	impervious	| NDVI	| elevation

**Target variable:**

ΔUHI_cell

**Predictor variables:**

Built environment characteristics associated with each grid cell.
