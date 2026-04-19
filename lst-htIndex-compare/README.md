# LST vs Heat Exposure Index Comparison Sub-Project

## Overview

This sub-project is part of the main **[Urban Heat Island Response to Extreme Heat Events](https://github.com/cchen744/uhi-extreme-heat-response/)** research initiative. It focuses on understanding how different thermal indicators (Land Surface Temperature and Heat Exposure Index) capture urban thermal variations from the perspective of **human living environment**.

## Motivation

The main project investigates UHI intensity during extreme heat events across multiple cities. However, different thermal metrics capture urban heat at different scales and perspectives:

- **Land Surface Temperature (LST)**: Satellite-based measure of Earth's surface thermal radiation; primarily reflects material properties (concrete, asphalt, vegetation)
- **Heat Exposure Index (HTIndex)**: Human-centric metric combining temperature, humidity, and apparent temperature; directly relates to human thermal comfort and physiological stress

This sub-project addresses the question: **Are LST and HTIndex spatially consistent?** Or do they reveal fundamentally different patterns of urban heat exposure relevant to human well-being?

## Key Research Questions

1. **Spatial Distribution**: Do LST and HTIndex hotspots align geographically?
2. **Magnitudes**: Which metric shows stronger spatial clustering (Moran's I)?
3. **Temporal Coherence**: Do seasonal and annual trends match across metrics?
4. **Climate Zone Differences**: How do humid (Miami) vs. arid (Phoenix) cities differ in LST-HTIndex relationships?
5. **Policy Implications**: Which metric better represents human thermal exposure for urban planning?

## Data & Methodology

### Current Phase: LST Analysis (Daytime, 2016-2026)

**Data Source**: 
- MODIS MOD11A1 (Terra satellite, daytime LST)
- 1 km resolution
- Seasonal aggregation (DJF, MAM, JJA, SON)
- Quality-assured pixels only (QA filtering in Google Earth Engine)

**Spatial Analysis**:
- **Moran's I**: Measures spatial autocorrelation (clustering of similar temperatures)
- **Hotspot Detection**: Identifies high-temperature areas (>75th percentile)
- **Descriptive Statistics**: Mean, median, std dev, spatial variability

**Study Cities**:
- **Miami** (FL): Humid subtropical, coastal, moderate LST range
- **Phoenix** (AZ): Hot desert, inland, extreme LST range
- Framework allows parametrized analysis for all 7 cities in the main project

### Next Phase: Heat Exposure Index

Once LST pipeline is validated, will add:
- Apparent temperature calculations
- Humidity data integration
- Heat stress indices (WBGT, Humidex, etc.)
- Cross-metric spatial correlation analysis

## Repository Structure

```
lst-htIndex-compare/
├── README.md                    # This file
├── lst_pipeline.py              # Production pipeline (GEE + analysis)
├── LST_analysis.ipynb           # Demo notebook (Miami vs Phoenix)
├── outputs/
│   ├── Miami/
│   │   ├── Miami_statistics.csv
│   │   └── Miami_mean_LST.png
│   ├── Phoenix/
│   │   ├── Phoenix_statistics.csv
│   │   └── Phoenix_mean_LST.png
│   ├── LST_combined_statistics.csv
│   ├── temperature_trends_comparison.png
│   ├── morans_i_comparison.png
│   └── [other visualizations]
```

## Usage

### Quick Start

```python
from lst_pipeline import create_config, LSSTPipeline

# Run for any city with parametrized config
config = create_config('Miami', start_year=2016, end_year=2026)
pipeline = LSSTPipeline(config)
results = pipeline.run()

# Results include statistics and visualizations
```

### Colab Setup

```bash
# Install dependencies
!pip install ee geemap pandas numpy xarray rasterio matplotlib seaborn scipy esda

# Authenticate GEE (first run only)
import ee
ee.Authenticate()
ee.Initialize()
```

## Key Findings (LST Phase)

**Temperature Differences** (2016-2026 mean):
- Miami: ~25.5°C
- Phoenix: ~35.8°C
- Difference: ~10.3°C (desert vs. coastal climate effect)

**Spatial Clustering (Moran's I)**:
- Phoenix shows stronger spatial autocorrelation → more clustered temperature patterns
- Miami shows weaker clustering → more heterogeneous local conditions

**Variability**:
- Phoenix: Higher spatial std dev (~8-10°C within city)
- Miami: Lower spatial std dev (~3-5°C within city)
- Suggests Phoenix has more pronounced UHI contrasts

**Hotspot Areas** (>75th percentile):
- Both cities show seasonal variation in high-temperature zones
- Phoenix summer (JJA) shows more extensive hotspots than Miami

## Connection to Main Project

This sub-project **validates LST as a proxy for human thermal exposure** within the main UHI study:

1. **Scale Alignment**: Confirms MODIS 1 km resolution sufficient for intra-urban analysis
2. **Temporal Consistency**: LST seasonal patterns align with known heat-health vulnerabilities
3. **Comparative Framework**: Establishes baseline for comparing LST with HTIndex
4. **Cities Coverage**: Miami and Phoenix are subset of 7-city main project cohort

Results will inform which metrics best represent urban thermal exposure in the main project's analysis of extreme heat events and UHI intensity.

## Human Living Environment Perspective

Rather than treating LST as abstract "surface temperature," this sub-project explicitly frames analysis in terms of **human thermal exposure**:

- **Daytime LST**: Reflects daytime outdoor conditions affecting pedestrians, outdoor workers
- **Spatial Clustering**: Identifies neighborhoods with persistent heat stress
- **Hotspots**: Reveals where populations experience highest thermal loads
- **City Comparisons**: Shows climate-dependent patterns relevant to adaptation strategies

Future HTIndex work will directly assess **apparent temperature** (how hot humans perceive it), making the human component explicit.

## Publications & References

Main project publications/preprints will be linked here as available.

## Contributing

This branch is under active development. For questions, issues, or contributions to the LST-HTIndex comparison framework, contact the project lead.

---

**Sub-project Lead**: Albert Chen (cchen744@wisc.edu)  
**Main Project**: Urban Heat Island Response to Extreme Heat Events  
**Last Updated**: April 2026
