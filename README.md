# Urban Heat Island Response to Extreme Heat Events

## 1. Research Question
a. How does urban heat island (UHI) intensity respond under extreme heat conditions compared to baseline summer conditions?
<br>b. Is the condition-dependence of UHI associated with differences in urban land use composition or other built-environment factors across cities?

## 2. Study Area and Data
- Daily PRISM temperature data (tmean, tmax)
- Urban area boundaries (Census urban areas)
- Study area: Chicago metropolitan region

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
