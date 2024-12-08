import geopandas as gpd
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
from pathlib import Path
import logging

# Set up logging
log_path = Path('../logs/step2_debugging.log')
log_path.parent.mkdir(parents=True, exist_ok=True)
logging.basicConfig(filename=log_path, level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Assuming occurence and insects are already loaded as GeoDataFrames
occurrence = pd.read_csv('/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_co-occurrences_1986-2012.csv')
out_path = '/home/goldma34/fire_insect_co-occurence/data/outputs/on/'



def process_fire_ids(fire_ids_chunk):
    results = []
    
    for fire_id in fire_ids_chunk:
        logging.info(f"fire id: {fire_id}")

        # Subset the occurence GeoDataFrame to the current Fire_ID
        occurrence_subset = occurrence[occurrence['Fire_ID'] == fire_id]

        # calculate the time since defoliation
        time_since_defoliation = (occurrence_subset['Fire_Year'] - occurrence_subset['Insect_Year'].max()).max()

        # calculate the number of insect years for each fire, 
        # but if TSD is 0 don't count the year closest to the fire year
        if time_since_defoliation == 0:
            cumulative_years = (occurrence_subset['Insect_Year'].nunique() - 1)
        else:
            cumulative_years = (occurrence_subset['Insect_Year'].nunique())

        # Calculate the total area of the fire that was overlapped.
        max_area = occurrence_subset['Overlap_Area'].max()
        max_percent = occurrence_subset['Percent_Overlap'].max()

        # append results
        results.append({'Fire_ID': fire_id, 
                        'Fire_Year': occurrence_subset['Fire_Year'].max(),
                        'Time_Since_Defoliation': time_since_defoliation, 
                        'Cumulative_Years': cumulative_years, 
                        'Max_Overlap_Area': max_area,
                        'Max_Overlap_Percent': max_percent})
    return results

if __name__ == '__main__':
    # Split the unique Fire_ID values into 20 chunks
    fire_ids = occurrence['Fire_ID'].unique()
    n_chunks = 20
    fire_id_chunks = np.array_split(fire_ids, n_chunks)

    # Use joblib to process the chunks in parallel
    results = Parallel(n_jobs=5)(delayed(process_fire_ids)(chunk) for chunk in fire_id_chunks)

    # Flatten the list of results
    results_flat = [item for sublist in results for item in sublist]

    # Convert the results to a DataFrame
    final_df = pd.DataFrame(results_flat)

    # Save the results to a CSV file
    final_df.to_csv(Path(out_path) / 'on_co-occurrence_defoliation_history.csv', index=False)