import pandas as pd
import random
import tempfile

# Read in fire-insect co-occurrence data
occurrence = pd.read_csv('/home/goldma34/fire_insect_co-occurence/data/outputs/qc/qc_fire_insect_co-occurrence_v1.csv')

# Select a random fire
random_fire = occurrence.sample(n=1)

# Create a temporary file
with tempfile.NamedTemporaryFile(delete=False, suffix='.csv') as temp_file:
    temp_file_path = temp_file.name
    # Save the selected fire to the temporary file
    random_fire.to_csv(temp_file_path, index=False)

# Print the path to the temporary file
print(temp_file_path)

