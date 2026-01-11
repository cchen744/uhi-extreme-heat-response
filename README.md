# uhi-extreme-heat-response
# Urban Heat Island Response to Extreme Heat Events

## 1. Research Question
a. How does urban heat island (UHI) intensity respond under extreme heat conditions compared to baseline summer conditions?
b. Is the condition-dependence of UHI associated with differences in urban land use composition or other built-environment factors across cities?

## 2. Study Area and Data
- Daily PRISM temperature data (tmean, tmax)
- Urban area boundaries (Census urban areas)
- Study area: Chicago metropolitan region

## 3. Method Overview
- Pixel-level urban vs non-urban temperature contrasts
- Apparent temperature–based extreme heat event definition
- Event-based comparison of UHI intensity

## 4. Repository Structure
uhi-extreme-heat-response/
│
├── README.md
├── environment.yml
│
├── data/
│   ├── raw/                # untouched source data 
│   ├── interim/            # masked, clipped, joined outputs
│   └── processed/          # analysis-ready tables/rasters
│
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_masking_and_aggregation.ipynb
│   ├── 03_uhi_construction.ipynb
│   └── 04_extreme_heat_events.ipynb
│
├── src/
│   ├── config.py
│   ├── masking.py
│   ├── temperature_metrics.py
│   ├── uhi.py
│   └── extreme_events.py
│
├── results/
│   ├── figures/
│   └── tables/
│
└── docs/
    └── methodology.md


## 5. Reproducibility
