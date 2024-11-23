import os
import sys
import logging
import concurrent.futures
import geopandas as gpd
import multiprocessing
from functools import partial
from shapely.validation import make_valid
from pathlib import Path

# Add path
sys.path.append(os.path.abspath("/home/goldma34/ch3-null-model/src/"))
from paths import *

# Configure logging for intersections
intersection_log_path = os.path.join("/home/goldma34/ch3-null-model/logs/", "intersection_log.log")
logging.basicConfig(filename=intersection_log_path, level=logging.INFO, format='%(asctime)s %(message)s')
error_handler = logging.FileHandler('/home/goldma34/ch3-null-model/logs/error.log')
error_handler.setLevel(logging.ERROR)
error_handler.setFormatter(logging.Formatter('%(asctime)s %(levelname)s:%(message)s'))
logging.getLogger().addHandler(error_handler)

def validate_and_fix_geometry(geometry):
    if geometry is None:
        return None
    if not geometry.is_valid:
        geometry = geometry.buffer(0)
    return geometry

def process_fire(fire, insects):
    try:
        # Create an empty list to store the results
        results = []

        # Fix geometry
        fire.geometry = validate_and_fix_geometry(fire.geometry)
        
        # Check if fire geometry is valid
        if fire.geometry is None:
            logging.error(f"Invalid geometry for fire_id {fire['Fire_ID']}")
            return results

        # Loop through each row in the insects GeoDataFrame
        for i, insects_row in insects.iterrows():
            # Fix geometry
            insects_row.geometry = validate_and_fix_geometry(insects_row.geometry)

            # Check if insects geometry is valid
            if insects_row.geometry is None:
                logging.error(f"Invalid geometry for insect row {i}")
                continue

            # Check if the geometries intersect
            if insects_row.geometry.intersects(fire.geometry):
                # Clip the intersection
                intersection = insects_row.geometry.intersection(fire.geometry)
                # Append the result to the results list
                results.append({
                    'year': insects_row['Year'],
                    'fire_id': fire['Fire_ID'],
                    'fire_year': fire['Fire_Year'],
                    'defoliated': 'defoliated',
                    'geometry': intersection
                })
                # Log the intersection
                logging.info(f"Intersection happened for fire_id {fire['Fire_ID']}")
            else:
                # Handle the case where there is no intersection
                logging.info(f"No intersection for fire_id {fire['Fire_ID']}")
        
        logging.info(f"Results for fire_id {fire['Fire_ID']}: {results}")
        return results
    except Exception as e:
        logging.error(f"Error processing fire {fire['Fire_ID']}: {e}")
        return []

def process_fires_parallel(fires, insects):
    all_results = []

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_fire, fire, insects) for _, fire in fires.iterrows()]
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                logging.info(f"Future result: {result}")
                all_results.extend(result)
            except Exception as e:
                logging.error(f"Error in future result: {e}")

    if all_results:
        return gpd.GeoDataFrame(all_results)
    else:
        return None

# Example usage

fire_path = fire_in_path
insect_path = insect_in_path
out_path = shp_out_path

# Load the fire shapefile
fires = gpd.read_file(fire_path)

# Convert 'Fire_Year' to integer
fires['Fire_Year'] = fires['Fire_Year'].astype(int)

# Ensure it's a GeoDataFrame and validate geometries
fires = gpd.GeoDataFrame(fires)
fires.loc[:, 'geometry'] = fires['geometry'].apply(lambda geom: validate_and_fix_geometry(geom))

# Filter fires by year
fires_filtered = fires[(fires['Fire_Year'] >= 1985) & (fires['Fire_Year'] <= 2012)]

# Fires
fires_filtered = gpd.GeoDataFrame(fires_filtered)

# Randomly select 30 fires
# TODO: Remove this line later, it's temporary for testing purposes
if len(fires_filtered) > 30:
    fires_filtered = fires_filtered.sample(n=30, random_state=1)

# Load the insects shapefile and ensure the same CRS
insects = gpd.read_file(insect_path)
insects = insects.to_crs(fires.crs)

result_gdf = process_fires_parallel(fires_filtered, insects)
final_shp_out_path = Path(out_path)
overlap_out_path = str(final_shp_out_path / "qc_fire_overlap_cases.shp")

if result_gdf is not None and not result_gdf.empty:
    result_gdf.to_file(overlap_out_path)
else:
    logging.info("No intersections found.")