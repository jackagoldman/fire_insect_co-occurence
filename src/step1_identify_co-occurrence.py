import os
import sys
import logging
import geopandas as gpd
from multiprocessing import Pool
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import pandas as pd


# Add path
sys.path.append(os.path.abspath("/home/goldma34/fire_insect_co-occurence/src/"))
from paths import *

#filepaths
fire_path = fire_in_path
insect_path = insect_in_path
out_path = shp_out_path

# Set up logging
log_path = Path('logs/step1_debugging.log')
log_path.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=log_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Define validate geometry function
def validate_and_fix_geometry(geometry):
    if not geometry.is_valid:
        geometry = geometry.buffer(0)
    return geometry if geometry.is_valid else None

########## indentify overlap cases for each fire
def fire_overlap_cases(fire, insects):
    # Empty results list
    results = []

    fire_year = int(fire['Fire_Year'].values[0])
    fire_id = fire['Fire_ID'].values[0]

    # check if column is year
    if 'Year' not in insects.columns:
    # Check if there is a column named 'year' and rename it to 'Year'
        if 'year' in insects.columns:
            insects = insects.rename(columns={'year': 'Year'})
            print("Column 'year' has been renamed to 'Year'.")
        else:
            print("Neither 'Year' nor 'year' column found.")
    else:
        print("Column 'Year' already exists.")

    
    # Convert Year column to int if it's not already
    if insects['Year'].dtype != 'int64':
        insects['Year'] = insects['Year'].astype(int)
        print("Column 'Year' has been converted to int.")
   
    # Filter insect year
    sbw = insects[insects['Year'] <= fire_year]

    # Match CRS
    fire.set_crs(sbw.crs, inplace=True, allow_override=True)
    fire = fire.to_crs(epsg=32633)
    
    sbw = sbw.to_crs(fire.crs)

    print(f"Processing fire {fire_id} from year {fire_year}")

    for row in sbw.itertuples():
        sbw_year = row.Year
        sbw_geom = row.geometry

        # Create a GeoDataFrame for the current sbw row
        sbw_gdf = gpd.GeoDataFrame([{'geometry': sbw_geom, 'Year': sbw_year}], crs=sbw.crs)

        # Calculate intersection
        intersection = gpd.overlay(fire, sbw_gdf, how='intersection')

        if intersection.empty:
            continue

        # Calculate area of intersection
        intersection_area = intersection['geometry'].area.sum()
        # Calculate area of fire
        fire_area = fire['geometry'].area.sum()
        # Calculate percent overlap
        percent_overlap = (intersection_area / fire_area) * 100

        # Store results only if percent overlap is greater than 5%
        if percent_overlap > 5:
            # Create DataFrame for the current result
            df = pd.DataFrame({
                'Fire_ID': [fire_id],
                'Fire_Year': [fire_year],
                'Insect_Year': [sbw_year],
                'Fire_Area': [fire_area],
                'Overlap_Area': [intersection_area],
                'Percent_Overlap': [percent_overlap]
            })
            # Append the DataFrame to the results list
            results.append(df)

    # Concatenate all DataFrames in the results list if not empty
    if results:
        final_df = pd.concat(results, ignore_index=True)
        return final_df
    else:
        return pd.DataFrame()  # Return an empty DataFrame if no results


######### interate over fires
def iterate_fires(fire, insect):
    results =[]

    fire = fire[(fire['Fire_Year'] >= 1985) & (fire['Fire_Year'] <= 2012)]

    # Ensure 'Fire_Year' column exists before filtering
    if 'Fire_Year' not in fire.columns:
        raise KeyError("Column 'Fire_Year' not found in fire2 GeoDataFrame")
    

    fire = fire.to_crs(epsg=32633)


    # Change the CRS of gdf to match the CRS of insect
    insect = insect.to_crs(fire.crs)


    for i, fire_row in fire.iterrows():
        fire_gdf = gpd.GeoDataFrame(fire_row).T
        fire_gdf = fire_gdf.set_geometry('geometry')  # Ensure geometry is set
        result_df = fire_overlap_cases(fire_gdf, insect)
        if not result_df.empty:
            results.append(result_df)

    # Concatenate all results
    if results:
        final_df = pd.concat(results, ignore_index=True)
    else:
        final_df = pd.DataFrame()
    
    return final_df

################ Split function
def split_dataframe(df, n):
    return [df.iloc[i::n] for i in range(n)]


if __name__ == '__main__':

# read in files
    fire = gpd.read_file(fire_path)
    insect = gpd.read_file(insect_path)

# Split the fire GeoDataFrame into 4 parts
    fire_parts = split_dataframe(fire, 4)

    # Run iterate_fires on each part in parallel
    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(iterate_fires, part, insect) for part in fire_parts]
        results = [future.result() for future in futures]

    # Combine the results
    final_df = pd.concat(results, ignore_index=True)

    # Save the GeoDataFrame to a file
    final_shp_out_path = Path(out_path)
    overlap_out_path = str(final_shp_out_path / "on_co-occurrences_1986-2012.csv")
    final_df.to_csv(overlap_out_path, index=False)