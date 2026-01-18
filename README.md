# Urban Heat Island Response to Extreme Heat Events

## 1. Research Question
a. How does urban heat island (UHI) intensity respond under extreme heat conditions compared to baseline summer conditions?
<br>b. Is the condition-dependence of UHI associated with differences in urban land use composition or other built-environment factors across cities?

## 2. Study Area and Data
- Daily PRISM temperature data (tmean, tmax)
- Urban area boundaries (Census urban areas)
- Study area
<br>**Core Study Cities**

| City         | Climate & Heat Exposure Context                                | UHI Risk & Inequality Signal                  | Role in Study Design                  |
| ------------ | -------------------------------------------------------------- | --------------------------------------------- | ------------------------------------- |
| Phoenix      | Hot desert climate with frequent extreme heat events           | High UHI index; extensive NOAA UHI mapping    | Extreme-heat-dominated reference case |
| Houston      | Humid subtropical climate; high nighttime temperatures         | Strong UHI amplification and social disparity | Humid heat amplification contrast     |
| Los Angeles  | Mediterranean coastal climate; dry summers                     | Pronounced intra-urban heat inequality        | Dry-heat and coastal mechanism case   |
| Chicago      | Humid continental climate; low-frequency but severe heat waves | Large UHI disparities in extreme heat         | Condition-dependence test case        |
| Philadelphia | Humid continental, legacy urban form                           | Stable redlining-associated UHI signal        | Compact, older urban morphology case  |
| Nashville    | Transitional subtropical climate; rapid urban growth           | NOAA-mapped UHI heterogeneity                 | Rapid-expansion urban form case       |

**Optional**
| City          | Intended Use                                            |
| ------------- | ------------------------------------------------------- |
| Columbus      | Midwest intra-regional comparison                       |
| San Francisco | Coastal city with weaker extreme heat; negative control |
| Jacksonville  | Humid subtropical, low-density urban form               |

Cities were selected to ensure climatic, urban form, and observed UHI variation while guaranteeing multi-source data comparability. Specifically:
(1) They span diverse extreme-heat climate regimes (e.g., Phoenix's hot desert, Houston's humid subtropical) to test for condition-dependent UHI amplification; (2) All are part of NOAA's Urban Heat Island Mapping Campaigns, providing independent empirical UHI data; (3) Each shows strong intra-urban LST disparities, particularly between historically redlined and non-redlined areas, linking UHI response to land use differences; (4) They represent contrasting development patterns (e.g., compact legacy cities, rapidly expanding metros), enabling analysis of how built-environment structure mediates UHI under extreme heat.

<p>
<img src=https://upload.wikimedia.org/wikipedia/commons/7/77/K%C3%B6ppen_Climate_Types_US_50.png align="middle" width=600 />
</p>

## 3. Method Overview
- Pixel-level urban vs non-urban temperature contrasts
- Apparent temperature–based extreme heat event definition
- Event-based comparison of UHI intensity

## 4. Repository Structure
<pre>
uhi-extreme-heat-response/
<br>├── README.md
<br>├── environment.yml
<br>├── data/
<br>   ├── raw/                # untouched source data 
<br>   ├── cities/            # masked, clipped, joined outputs
<br>   └── processed/          # analysis-ready tables/rasters
<br>├── notebooks/
<br>  ├── 01_data_exploration.ipynb
<br>   ├── 02_uhi_extreme_heat_response.ipynb
<br>   ├── 03_uhi_construction.ipynb
<br>  └── 04_extreme_heat_events.ipynb
<br>├── src/
<br>  ├── config.py
<br>  ├── masking.py
<br>  ├── temperature_metrics.py
<br>  ├── uhi.py
<br>  └── extreme_events.py
<br>├── results/
<br>   ├── figures/
<br>   └── tables/
<br>└── docs/
<br>    └── methodology.md
uhi_pipeline.py              # preprocessing pipeline that extracts and cleans raw data using google earth engine
</pre>

## 5. Reproducibility
<text>
