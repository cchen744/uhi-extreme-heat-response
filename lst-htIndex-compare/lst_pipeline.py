"""
LST Summer (JJA) Analysis Pipeline
Focused on June-July-August spatial temperature distribution
MODIS MOD11A1 daytime LST, 2020-2024
"""

import ee
import pandas as pd
import numpy as np
import os
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# Optional: for map visualization
try:
    import folium
    import rasterio
    HAS_VIZ = True
except ImportError:
    HAS_VIZ = False


# ============================================================================
# CONFIG CLASS
# ============================================================================

class LSTConfig:
    """Configuration for JJA LST analysis"""

    def __init__(self, city_name, city_bounds, start_year=2020, end_year=2024, output_dir='./outputs'):
        """
        Args:
            city_name (str): Name of city
            city_bounds (list): [min_lon, min_lat, max_lon, max_lat]
            start_year (int): Start year (default 2020)
            end_year (int): End year (default 2024, inclusive)
            output_dir (str): Output directory path
        """
        self.city_name = city_name
        self.city_bounds = city_bounds
        self.start_year = start_year
        self.end_year = end_year
        self.output_dir = Path(output_dir) / city_name
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # MODIS settings
        self.modis_product = 'MODIS/061/MOD11A1'  # Terra daytime LST
        self.lst_band = 'LST_Day_1km'
        self.qa_band = 'QC_Day'

    def create_geom(self):
        """Create GEE geometry from bounds"""
        min_lon, min_lat, max_lon, max_lat = self.city_bounds
        return ee.Geometry.Rectangle([min_lon, min_lat, max_lon, max_lat])


# ============================================================================
# GEE DATA FETCHER
# ============================================================================

class GEELSTFetcher:
    """Fetch MODIS LST data from GEE"""

    def __init__(self):
        try:
            ee.Initialize()
        except Exception as e:
            print(f"GEE initialization note: {e}")

    def fetch_jja_data(self, config):
        """
        Fetch JJA (June-July-August) LST data for each year
        Returns dict: {year: ee.Image}
        """
        jja_data = {}
        geom = config.create_geom()

        for year in range(config.start_year, config.end_year + 1):
            # JJA: months 6, 7, 8
            start_date = ee.Date.fromYMD(year, 6, 1)  # June 1
            end_date = ee.Date.fromYMD(year, 9, 1)    # Sept 1 (exclusive)

            collection = (ee.ImageCollection(config.modis_product)
                          .filterDate(start_date, end_date)
                          .filterBounds(geom)
                          .select([config.lst_band, config.qa_band]))

            # Apply QA filtering
            filtered = collection.map(lambda img: self._qa_filter(img, config))

            # Compute median across 3 months
            median_img = filtered.select(config.lst_band).median()

            jja_data[year] = median_img
            print(f"  {year} JJA: fetched")

        return jja_data

    @staticmethod
    def _qa_filter(image, config):
        """QA filtering for MOD11A1: keep quality 0 and 1"""
        qa = image.select(config.qa_band)
        good_quality = qa.lte(1)  # 0=good, 1=marginal
        return image.select(config.lst_band).updateMask(good_quality)


# ============================================================================
# SPATIAL ANALYSIS
# ============================================================================

class SpatialAnalyzer:
    """Compute spatial statistics"""

    @staticmethod
    def compute_stats(image, region, scale=1000):
        """Compute basic stats on GEE image"""
        stats = image.reduceRegion(
            reducer=ee.Reducer.mean().combine(
                ee.Reducer.median(),
                sharedInputs=True
            ).combine(
                ee.Reducer.stdDev(),
                sharedInputs=True
            ).combine(
                ee.Reducer.minMax(),
                sharedInputs=True
            ),
            geometry=region,
            scale=scale
        ).getInfo()

        return stats

    @staticmethod
    def identify_hotspots(data_array, threshold_percentile=75):
        """Identify hotspot pixels (high temperature)"""
        threshold = np.nanpercentile(data_array, threshold_percentile)
        hotspots = data_array > threshold
        return hotspots, threshold


# ============================================================================
# VISUALIZATION & EXPORT
# ============================================================================

class LSTVisualizer:
    """Create visualizations and maps"""

    @staticmethod
    def export_to_geotiff(image, config, scale=500):
        """Export GEE image to GeoTIFF in Google Drive"""
        geom = config.create_geom()

        # Unscale to Celsius
        lst_celsius = image.multiply(0.02).subtract(273.15)

        task = ee.batch.Export.image.toDrive(
            image=lst_celsius,
            description=f'{config.city_name}_JJA_LST_{config.start_year}_{config.end_year}',
            folder='LST_maps',
            scale=scale,
            region=geom,
            fileFormat='GeoTIFF'
        )
        task.start()
        print(f"Exporting {config.city_name} JJA LST to Google Drive...")
        return task

    @staticmethod
    def create_hotspot_map(geotiff_path, city_name, output_html, city_center):
        """
        Create interactive Folium map with hotspots

        Args:
            geotiff_path: Path to downloaded GeoTIFF
            city_name: Name of city
            output_html: Output HTML path
            city_center: [lat, lon]
        """
        if not HAS_VIZ:
            print("Install rasterio and folium for map visualization")
            return

        # Read GeoTIFF
        with rasterio.open(geotiff_path) as src:
            data = src.read(1)
            bounds = src.bounds

        # Identify hotspots
        threshold = np.nanpercentile(data, 75)
        hotspot_mask = (data > threshold).astype(float)

        # Normalize data for coloring
        norm_lst = Normalize(vmin=np.nanpercentile(data, 2),
                             vmax=np.nanpercentile(data, 98))
        norm_hotspot = Normalize(vmin=0, vmax=1)

        cmap_lst = cm.get_cmap('RdYlBu_r')
        cmap_hotspot = cm.get_cmap('Reds')

        # Convert to RGB images
        lst_rgb = cmap_lst(norm_lst(data))
        hotspot_rgb = cmap_hotspot(norm_hotspot(hotspot_mask))

        # Save as PNG
        import tempfile
        with tempfile.TemporaryDirectory() as tmpdir:
            lst_png = f'{tmpdir}/lst.png'
            hotspot_png = f'{tmpdir}/hotspot.png'

            plt.imsave(lst_png, lst_rgb)
            plt.imsave(hotspot_png, hotspot_rgb)

            # Create Folium map
            m = folium.Map(
                location=city_center,
                zoom_start=11,
                tiles='CartoDB positron'
            )

            # Add LST layer
            folium.raster_layers.ImageOverlay(
                image=lst_png,
                bounds=[[bounds.bottom, bounds.left], [bounds.top, bounds.right]],
                opacity=0.7,
                name='LST (°C)',
                show=True
            ).add_to(m)

            # Add hotspot overlay
            folium.raster_layers.ImageOverlay(
                image=hotspot_png,
                bounds=[[bounds.bottom, bounds.left], [bounds.top, bounds.right]],
                opacity=0.4,
                name='Hotspots (>75th percentile)',
                show=False
            ).add_to(m)

            # Add layer control
            folium.LayerControl().add_to(m)

            m.save(output_html)
            print(f"Map saved to {output_html}")


# ============================================================================
# MAIN PIPELINE
# ============================================================================

class JJALSTPipeline:
    """Main JJA LST analysis pipeline"""

    def __init__(self, config):
        self.config = config
        self.fetcher = GEELSTFetcher()
        self.analyzer = SpatialAnalyzer()
        self.visualizer = LSTVisualizer()
        self.results = {}

    def run(self, export_geotiff=True):
        """
        Execute full JJA LST pipeline

        Args:
            export_geotiff (bool): Export to GeoTIFF for map visualization
        """
        print(f"\n{'='*60}")
        print(f"JJA LST Analysis Pipeline: {self.config.city_name}")
        print(f"Period: {self.config.start_year}-{self.config.end_year}")
        print(f"Months: June, July, August")
        print(f"{'='*60}\n")

        # Step 1: Fetch JJA data
        print("Step 1: Fetching MODIS JJA LST data...")
        jja_data_gee = self.fetcher.fetch_jja_data(self.config)
        self.results['jja_data_gee'] = jja_data_gee

        # Step 2: Compute statistics
        print("Step 2: Computing statistics...")
        stats_list = self._compute_all_stats(jja_data_gee)
        stats_df = pd.DataFrame(stats_list)
        self.results['statistics'] = stats_df

        # Step 3: Export GeoTIFF for visualization
        if export_geotiff:
            print("Step 3: Exporting GeoTIFF...")
            jja_mean = self._compute_jja_mean(jja_data_gee)
            task = self.visualizer.export_to_geotiff(jja_mean, self.config)
            self.results['export_task'] = task
            print("  → Download from Google Drive when ready")

        # Step 4: Export statistics
        print("Step 4: Exporting statistics...")
        self._export_results(stats_df)

        print("\n✓ Pipeline complete!\n")
        return self.results

    def _compute_jja_mean(self, jja_data_gee):
        """Compute mean across all years"""
        images = [img for img in jja_data_gee.values() if img is not None]
        if not images:
            return None
        return ee.ImageCollection(images).mean()

    def _compute_all_stats(self, jja_data_gee):
        """Compute stats for each year"""
        records = []
        geom = self.config.create_geom()

        for year in sorted(jja_data_gee.keys()):
            image = jja_data_gee[year]

            if image is not None:
                # Get statistics from GEE
                stats = self.analyzer.compute_stats(image, geom)

                if stats and f'{self.config.lst_band}_mean' in stats:
                    # Convert from Kelvin to Celsius
                    mean_k = stats[f'{self.config.lst_band}_mean'] * 0.02
                    median_k = stats[f'{self.config.lst_band}_median'] * 0.02
                    std_k = stats[f'{self.config.lst_band}_stdDev'] * 0.02
                    min_k = stats[f'{self.config.lst_band}_min'] * 0.02
                    max_k = stats[f'{self.config.lst_band}_max'] * 0.02

                    records.append({
                        'year': year,
                        'season': 'JJA',
                        'mean_temp_c': mean_k - 273.15,
                        'median_temp_c': median_k - 273.15,
                        'std_temp_c': std_k,
                        'min_temp_c': min_k - 273.15,
                        'max_temp_c': max_k - 273.15,
                        'temp_range_c': (max_k - min_k)
                    })

        return records

    def _export_results(self, stats_df):
        """Export statistics to CSV"""
        output_file = self.config.output_dir / f"{self.config.city_name}_JJA_statistics.csv"
        stats_df.to_csv(output_file, index=False)
        print(f"  → Statistics: {output_file}")
        print(f"  → Output dir: {self.config.output_dir}")


# ============================================================================
# CITY PRESETS
# ============================================================================

CITY_PRESETS = {
    'Miami': {
        'bounds': [-80.3, 25.7, -80.1, 25.9],
        'center': [25.8, -80.2]
    },
    'Phoenix': {
        'bounds': [-112.15, 33.35, -112.0, 33.5],
        'center': [33.425, -112.075]
    }
}


def create_config(city_name, start_year=2020, end_year=2024, output_dir='./outputs'):
    """Create config from city preset"""
    if city_name not in CITY_PRESETS:
        raise ValueError(f"City {city_name} not in presets. Available: {list(CITY_PRESETS.keys())}")

    bounds = CITY_PRESETS[city_name]['bounds']
    return LSTConfig(city_name, bounds, start_year, end_year, output_dir)


# ============================================================================
# USAGE EXAMPLE
# ============================================================================

if __name__ == "__main__":
    # Example: Run for Miami
    config = create_config('Miami', start_year=2020, end_year=2024)
    pipeline = JJALSTPipeline(config)
    results = pipeline.run(export_geotiff=True)

    # Print results
    print(results['statistics'])
