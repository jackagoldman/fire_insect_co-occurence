import os
import sys
import logging
import concurrent.futures
import geopandas as gpd
import multiprocessing
from functools import partial
from shapely.validation import make_valid
from pathlib import Path
import pandas as pd
from shapely.geometry import shape
from concurrent.futures import ProcessPoolExecutor, as_completed

# Add path
sys.path.append(os.path.abspath("/home/goldma34/fire_insect_co-occurrence/src/"))
from paths import *

# Ensure the log directory exists
log_dir = "/home/goldma34/fire_insect_co-occurrence/logs/"
os.makedirs(log_dir, exist_ok=True)

# Configure logging for errors
error_log_path = os.path.join(log_dir, "error.log")
logging.basicConfig(filename=error_log_path, level=logging.ERROR, format='%(asctime)s %(levelname)s:%(message)s')


def validate_and_fix_geometry(geometry):
    if not geometry.is_valid:
        geometry = geometry.buffer(0)
    return geometry if geometry.is_valid else None

def process_fire(fire, insects):
    results = []
    try:
        # Filter out rows where the insect year is greater than the fire year
        filtered_insects = insects[insects['Year'] <= fire['Fire_Year']]
        

        for i, insects_row in filtered_insects.iterrows():
            insect_geometry = validate_and_fix_geometry(insects_row.geometry)
            if insect_geometry is None:
                logging.error(f"Invalid geometry for insect_id {insects_row['Insect_ID']}")
                continue

            if fire.geometry.intersects(insect_geometry):
                results.append({'Fire_ID': fire['Fire_ID'], 
                                'Insect_Year': insects_row['Year'],
                                'Fire_Year': fire['Fire_Year'],
                                'defoliated': 'defoliated',
                                'geometry': fire.geometry})


    except Exception as e:
        logging.error(f"Error processing fire_id {fire['Fire_ID']}: {e}")
    
    return results

def main(fire_gdf, insects_gdf):
    results = []
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_fire, fire, insects_gdf) for _, fire in fire_gdf.iterrows()]
        for future in as_completed(futures):
            results.extend(future.result())
    return results

if __name__ == "__main__":
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

    

    # Load the insects shapefile and ensure the same CRS
    insects = gpd.read_file(insect_path)
    insects = insects.to_crs(fires.crs)

    # Process fires in parallel
    results = main(fires_filtered, insects)

    # Convert results to GeoDataFrame
    result_gdf = gpd.GeoDataFrame(results)

    # Save the GeoDataFrame to a file
    final_shp_out_path = Path(out_path)
    overlap_out_path = str(final_shp_out_path / "qc_co-occurences_1986-2012.shp")
    result_gdf.to_file(overlap_out_path)

