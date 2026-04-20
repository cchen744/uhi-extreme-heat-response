"""
LST Pipeline for geospatial analysis of Land Surface Temperature
Fetches MODIS MOD11A1 data, performs seasonal aggregation, spatial analysis, and visualization
"""

from logging import config

import ee
import geemap
import pandas as pd
import numpy as np
import xarray as xr
import rasterio
from rasterio.transform import from_bounds
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import Normalize
import seaborn as sns
from esda.moran import Moran
from scipy.spatial import distance_matrix
from scipy import stats
import warnings
import os
from datetime import datetime
from pathlib import Path

warnings.filterwarnings('ignore')

# ============================================================================
# CONFIG CLASS
# ============================================================================

class LSTConfig:
    """Configuration for LST analysis pipeline"""

    def __init__(self, city_name, city_bounds, start_year=2016, end_year=2026, output_dir='./outputs'):
        """
        Args:
            city_name (str): Name of city for output labeling
            city_bounds (list): [min_lon, min_lat, max_lon, max_lat]
            start_year (int): Start year
            end_year (int): End year (inclusive)
            output_dir (str): Output directory path
        """
        self.city_name = city_name
        self.city_bounds = city_bounds  # [min_lon, min_lat, max_lon, max_lat]
        self.start_year = start_year
        self.end_year = end_year
        self.output_dir = Path(output_dir) / city_name
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # MODIS settings
        self.modis_product = 'MODIS/061/MOD11A1'  # Terra daytime LST
        self.lst_band = 'LST_Day_1km'
        self.qa_band = 'QC_Day'

        # Seasonal definitions (month ranges)
        self.seasons = {
            'DJF': [12, 1, 2],    # Winter
            'MAM': [3, 4, 5],     # Spring
            'JJA': [6, 7, 8],     # Summer
            'SON': [9, 10, 11]    # Fall
        }

    def create_geom(self):
        """Create GEE geometry from bounds or TIGER"""
        if hasattr(self, '_tiger_geom'):
            return self._tiger_geom
        else:
            min_lon, min_lat, max_lon, max_lat = self.city_bounds
            return ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])


# ============================================================================
# GEE DATA FETCHER
# ============================================================================

class GEELSTFetcher:
    """Fetch and process MODIS LST data from Google Earth Engine"""

    def __init__(self):
        try:
            ee.Initialize()
        except:
            print("GEE not initialized. Run ee.Authenticate() first if needed.")

    def fetch_seasonal_data(self, config):
        """
        Fetch MODIS LST data for each season and year
        Returns dict: {year: {season: ee.Image}}
        """
        seasonal_data = {}
        geom = config.create_geom()

        for year in range(config.start_year, config.end_year + 1):
            seasonal_data[year] = {}

            for season_name, months in config.seasons.items():
                images = []

                for month in months:
                    # Handle year boundary for DJF (December of previous year)
                    if month == 12:
                        y = year - 1
                    else:
                        y = year

                    start_date = ee.Date.fromYMD(y, month, 1)
                    end_date = start_date.advance(1, 'month')

                    # Filter MODIS data
                    collection = (ee.ImageCollection(config.modis_product)
                                  .filterDate(start_date, end_date)
                                  .filterBounds(geom)
                                  .select([config.lst_band, config.qa_band]))

                    # Apply QA filtering (keep only good quality pixels)
                    filtered = collection.map(lambda img: self._qa_filter(img, config))
                    images.append(filtered)

                # Merge monthly images and compute median
                merged = ee.ImageCollection(images).flatten()

                if merged.size().getInfo() > 0:
                    median_img = merged.select(config.lst_band).median()
                    seasonal_data[year][season_name] = median_img

        return seasonal_data

    @staticmethod
    def _qa_filter(image, config):
        """
        QA filtering for MOD11A1
        QC_Day values: 0=good, 1=marginal, 2=snow/ice, 3=cloud
        Keep only good quality (0) and marginal (1)
        """
        qa = image.select(config.qa_band)
        good_quality = qa.lte(1)  # 0 or 1 are acceptable
        return image.select(config.lst_band).updateMask(good_quality)

    def export_to_geotiff(self, image, config, season, year, scale=1000):
        """
        Export GEE image to GeoTIFF
        LST is in Kelvin with scale factor of 0.02
        """
        geom = config.create_geom()
        filename = f"{config.city_name}_{year}_{season}_LST.tif"
        filepath = config.output_dir / filename

        # Unscale LST (convert from raw to Kelvin)
        lst_kelvin = image.multiply(0.02)

        # Convert to Celsius
        lst_celsius = lst_kelvin.subtract(273.15)

        # Export
        task = ee.batch.Export.image.toDrive(
            image=lst_celsius,
            description=f"LST_{config.city_name}_{year}_{season}",
            region=geom,
            scale=scale,
            fileFormat='GeoTIFF'
        )
        task.start()
        return task


# ============================================================================
# SPATIAL ANALYSIS
# ============================================================================

class SpatialAnalyzer:
    """Compute spatial statistics and indices"""

    @staticmethod
    def compute_moran_i(data_array, mask=None):
        """
        Compute Moran's I for spatial autocorrelation
        Higher values = stronger spatial clustering

        Args:
            data_array (np.array): 2D array of values
            mask (np.array): Boolean mask (True = valid data)

        Returns:
            dict: {'morans_i': float, 'p_value': float}
        """
        # Flatten and handle masked values
        if mask is not None:
            valid_data = data_array[mask]
            valid_coords = np.argwhere(mask)
        else:
            valid_data = data_array.flatten()
            valid_coords = np.argwhere(~np.isnan(data_array))

        if len(valid_data) < 3:
            return {'morans_i': np.nan, 'p_value': np.nan}

        # Compute spatial weights matrix (inverse distance)
        distances = distance_matrix(valid_coords, valid_coords)
        # Avoid division by zero
        np.fill_diagonal(distances, 1)
        weights = 1.0 / (distances + 1e-6)
        np.fill_diagonal(weights, 0)

        # Normalize weights
        weights = weights / weights.sum(axis=1, keepdims=True)

        # Compute Moran's I
        moran = Moran(valid_data, weights)

        return {
            'morans_i': moran.I,
            'p_value': moran.p_norm,
            'z_score': (moran.I - moran.EI) / np.sqrt(moran.VI_norm)
        }

    @staticmethod
    def compute_basic_stats(data_array, mask=None):
        """Compute basic statistics"""
        if mask is not None:
            data = data_array[mask]
        else:
            data = data_array[~np.isnan(data_array)]

        return {
            'mean': np.nanmean(data),
            'median': np.nanmedian(data),
            'std': np.nanstd(data),
            'min': np.nanmin(data),
            'max': np.nanmax(data),
            'count': len(data)
        }

    @staticmethod
    def detect_hotspots(data_array, threshold_percentile=75):
        """Identify hotspots (high temperature areas)"""
        threshold = np.nanpercentile(data_array, threshold_percentile)
        hotspots = data_array > threshold
        return hotspots, threshold


# ============================================================================
# VISUALIZATION
# ============================================================================

class LSTVisualizer:
    """Create visualizations for LST analysis"""

    @staticmethod
    def static_map(data_array, title, output_path, cmap='RdYlBu_r', vmin=None, vmax=None):
        """Create static LST map with mean/median values"""
        fig, ax = plt.subplots(figsize=(10, 8), dpi=150)

        if vmin is None:
            vmin = np.nanpercentile(data_array, 2)
        if vmax is None:
            vmax = np.nanpercentile(data_array, 98)

        im = ax.imshow(data_array, cmap=cmap, vmin=vmin, vmax=vmax,
                       origin='upper', interpolation='nearest')

        cbar = plt.colorbar(im, ax=ax, label='Temperature (°C)')
        ax.set_title(title, fontsize=14, fontweight='bold')
        ax.set_xlabel('Longitude (pixels)')
        ax.set_ylabel('Latitude (pixels)')

        plt.tight_layout()
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()

        return output_path

    @staticmethod
    def create_animation(data_dict, output_path, title_base, cmap='RdYlBu_r'):
        """
        Create animated time series of LST

        Args:
            data_dict (dict): {year: {season: np.array}}
            output_path (str): Output path for GIF/MP4
            title_base (str): Base title for frames
            cmap (str): Colormap
        """
        frames = []
        years = sorted(data_dict.keys())
        seasons = ['DJF', 'MAM', 'JJA', 'SON']

        # Collect all data for vmin/vmax normalization
        all_data = []
        for year in years:
            for season in seasons:
                if season in data_dict[year] and data_dict[year][season] is not None:
                    all_data.append(data_dict[year][season])

        if not all_data:
            print("No valid data for animation")
            return

        all_data = np.concatenate([d.flatten() for d in all_data])
        vmin = np.nanpercentile(all_data, 2)
        vmax = np.nanpercentile(all_data, 98)

        fig, ax = plt.subplots(figsize=(10, 8), dpi=100)

        frame_count = 0
        for year in years:
            for season in seasons:
                if season in data_dict[year] and data_dict[year][season] is not None:
                    data = data_dict[year][season]

                    ax.clear()
                    im = ax.imshow(data, cmap=cmap, vmin=vmin, vmax=vmax,
                                  origin='upper', interpolation='nearest')
                    ax.set_title(f'{title_base} | {year} {season}',
                               fontsize=12, fontweight='bold')
                    if frame_count == 0:
                        plt.colorbar(im, ax=ax, label='Temperature (°C)')

                    fig.canvas.draw()
                    frame_count += 1

        plt.tight_layout()
        plt.savefig(output_path, bbox_inches='tight')
        plt.close()

        print(f"Animation frame blueprint saved to {output_path}")


# ============================================================================
# PIPELINE ORCHESTRATOR
# ============================================================================

class LSSTPipeline:
    """Main pipeline orchestrator"""

    def __init__(self, config):
        self.config = config
        self.fetcher = GEELSTFetcher()
        self.analyzer = SpatialAnalyzer()
        self.visualizer = LSTVisualizer()
        self.results = {}

    def run(self, export_geotiff=False):
        """
        Execute full pipeline

        Args:
            export_geotiff (bool): Whether to export intermediate GeoTIFFs to Drive
        """
        print(f"\n{'='*60}")
        print(f"LST Pipeline: {self.config.city_name}")
        print(f"Period: {self.config.start_year}-{self.config.end_year}")
        print(f"{'='*60}\n")

        # Step 1: Fetch seasonal data
        print("Step 1: Fetching MODIS data...")
        seasonal_data_gee = self.fetcher.fetch_seasonal_data(self.config)

        # Step 2: Convert to local arrays and analyze
        print("Step 2: Processing and analyzing data...")
        seasonal_data_local = self._gee_to_local(seasonal_data_gee)
        self.results['seasonal_data'] = seasonal_data_local

        # Step 3: Compute statistics
        print("Step 3: Computing spatial statistics...")
        stats_df = self._compute_all_stats(seasonal_data_local)
        self.results['statistics'] = stats_df

        # Step 4: Create visualizations
        print("Step 4: Creating visualizations...")
        self._create_all_visualizations(seasonal_data_local)

        # Step 5: Export results
        print("Step 5: Exporting results...")
        self._export_results(stats_df)

        print("\n✓ Pipeline complete!\n")
        return self.results

    def _gee_to_local(self, seasonal_data_gee):
        """Convert GEE images to numpy arrays"""
        seasonal_data_local = {}
        geom = self.config.create_geom()

        for year, seasons_dict in seasonal_data_gee.items():
            seasonal_data_local[year] = {}

            for season, image in seasons_dict.items():
                try:
                    # Sample image to numpy array
                    sample = image.sample(scale=1000, numPixels=10000, geometries=True)
                    data_list = sample.aggregate_list('LST_Day_1km').getInfo()

                    if data_list:
                        # Convert to array
                        arr = np.array(data_list)
                        # Remove outliers (scale factor applied, should be ~200-350K)
                        arr = np.array([x for x in arr if x is not None and 200 < x < 350]) * 0.02 - 273.15
                        seasonal_data_local[year][season] = arr
                except:
                    seasonal_data_local[year][season] = None

        return seasonal_data_local

    def _compute_all_stats(self, seasonal_data_local):
        """Compute Moran's I and basic statistics for all seasons"""
        records = []

        for year in sorted(seasonal_data_local.keys()):
            for season in ['DJF', 'MAM', 'JJA', 'SON']:
                if season in seasonal_data_local[year]:
                    data = seasonal_data_local[year][season]

                    if data is not None and len(data) > 0:
                        # Basic stats
                        basic_stats = self.analyzer.compute_basic_stats(data)

                        # Reshape to rough grid for Moran's I
                        side = int(np.sqrt(len(data)))
                        data_2d = data[:side**2].reshape(side, side)

                        # Moran's I
                        moran_stats = self.analyzer.compute_moran_i(data_2d)

                        # Hotspots
                        hotspots, threshold = self.analyzer.detect_hotspots(data_2d)

                        records.append({
                            'year': year,
                            'season': season,
                            'mean_temp': basic_stats['mean'],
                            'median_temp': basic_stats['median'],
                            'std_temp': basic_stats['std'],
                            'min_temp': basic_stats['min'],
                            'max_temp': basic_stats['max'],
                            'n_pixels': basic_stats['count'],
                            'morans_i': moran_stats['morans_i'],
                            'morans_i_pvalue': moran_stats['p_value'],
                            'morans_i_zscore': moran_stats['z_score'],
                            'hotspot_pct': (hotspots.sum() / hotspots.size) * 100
                        })

        return pd.DataFrame(records)

    def _create_all_visualizations(self, seasonal_data_local):
        """Create maps and animations"""
        # Static map: mean across entire period
        all_data = []
        for year in sorted(seasonal_data_local.keys()):
            for season in ['DJF', 'MAM', 'JJA', 'SON']:
                if season in seasonal_data_local[year] and seasonal_data_local[year][season] is not None:
                    all_data.append(seasonal_data_local[year][season])

        if all_data:
            overall_mean = np.mean([np.mean(d) for d in all_data])
            side = int(np.sqrt(len(all_data[0])))
            mean_grid = all_data[0][:side**2].reshape(side, side)

            title = f"{self.config.city_name} - Mean LST ({self.config.start_year}-{self.config.end_year})"
            output = self.config.output_dir / f"{self.config.city_name}_mean_LST.png"
            self.visualizer.static_map(mean_grid, title, str(output))

    def _export_results(self, stats_df):
        """Export statistics to CSV"""
        output_file = self.config.output_dir / f"{self.config.city_name}_statistics.csv"
        stats_df.to_csv(output_file, index=False)
        print(f"  → Statistics: {output_file}")
        print(f"  → Output dir: {self.config.output_dir}")


# ============================================================================
# CITY PRESETS
# ============================================================================

CITY_PRESETS = {
    'Miami': {
        'bounds': [-80.3, 25.7, -80.1, 25.9],  # Miami, Florida
        'climate': 'Humid Subtropical'
    },
    'Phoenix': {
        'bounds': [-112.15, 33.35, -112.0, 33.5],  # Phoenix, Arizona
        'climate': 'Hot Desert'
    }
}


def create_config(city_name, fips_or_preset=None, start_year=2020, end_year=2024, output_dir='./outputs'):
    """
    Create config from city preset OR FIPS code
    
    Usage:
        create_config('Miami')  # Uses preset bounds
        create_config('Miami', '1204000', 2020, 2024)  # Uses FIPS code
    """
    config = LSTConfig(city_name, city_bounds=[0,0,0,0], start_year=start_year, end_year=end_year, output_dir=output_dir)
    # If FIPS code provided, fetch geometry from TIGER
    if fips_or_preset and isinstance(fips_or_preset, str) and len(fips_or_preset) == 7:
        # FIPS code: fetch from TIGER
        state_fips = fips_or_preset[:2]
        place_fips = fips_or_preset[2:]
        
        tiger = ee.FeatureCollection('TIGER/2020/PLACES').filter(
            ee.Filter.eq('STATEFP', state_fips)
        ).filter(
            ee.Filter.eq('PLACEFP', place_fips)
        )
        config._tiger_geom = tiger.geometry()
    else:
        # Use preset bounds
        if city_name not in CITY_PRESETS:
            raise ValueError(f"City {city_name} not in presets. Available: {list(CITY_PRESETS.keys())}")
        config.city_bounds = CITY_PRESETS[city_name]['bounds']
    
    return config


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

if __name__ == "__main__":
    # Example: Run for Miami
    config = create_config('Miami', start_year=2016, end_year=2026)
    pipeline = LSSTPipeline(config)
    results = pipeline.run()

    print(results['statistics'].head())
