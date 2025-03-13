# DEFAULTS

import glob
from pathlib import Path, PurePath



## Step 1: file paths
## file path ##
fire_in_path ="/home/goldma34/fire_insect_co-occurence/data/inputs/on/fire/ON_FirePerimeters_85to2020_v00.shp"
insect_in_path = "/home/goldma34/fire_insect_co-occurence/data/inputs/on/sbw/sbw_60_20_v2.shp"

# outpath
shp_out_path = "/home/goldma34/fire_insect_co-occurence/data/outputs/on/"

# log path
error_log_path = "/home/goldma34/fire_insect_co-occurence/logs/"

# Step 3: file paths
# GEE paths
ontario_fires = "users/jandrewgoldman/Ont_BurnSeverity_Trends/ON_FirePerimeters_85to2020_v0"
quebec_fires = "users/jandrewgoldman/qc-fire-perims-shield-2"
