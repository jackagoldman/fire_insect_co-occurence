#!/bin/bash

# Run the select_random_fire.py script and capture the path to the temporary file
TEMP_FILE=$(python3 /home/goldma34/fire_insect_co-occurence/test/select_random_fire.py)

# Export the path to the temporary file as an environment variable
export TEMP_FIRE_PATH="$TEMP_FILE"

# Backup the original step3_landtrendr.py script
cp /home/goldma34/fire_insect_co-occurence/src/step3_landtrendr.py /home/goldma34/fire_insect_co-occurence/src/step3_landtrendr_backup.py

# Modify the step3_landtrendr.py script to read from the temporary file
sed -i '10i# Get the path to the temporary file from the environment variable\ntemp_file_path = os.getenv("TEMP_FIRE_PATH")\n\n# Read in the selected fire\nselected_fire = pd.read_csv(temp_file_path)\n' /home/goldma34/fire_insect_co-occurence/src/step3_landtrendr.py

# Run the step3_landtrendr.py script
python3 /home/goldma34/fire_insect_co-occurence/src/step3_landtrendr.py

# Revert the changes to the step3_landtrendr.py script
mv /home/goldma34/fire_insect_co-occurence/src/step3_landtrendr_backup.py /home/goldma34/fire_insect_co-occurence/src/step3_landtrendr.py

# Clean up the temporary file
rm "$TEMP_FILE"