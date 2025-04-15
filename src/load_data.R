# Required Data - prep 

library(dplyr)
library(sf)

# read in history data
history <- read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_defoliation_history_wx.csv")

# shorten host_percentage to host_pct
history <- history %>% 
  rename("host_pct" = "host_percentage")

# read in centroids
centroids <- read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_fire_centroids.csv")

# read in recovery
recovery_defol = read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_recovery_magnitude/on_recovery_magnitude.csv")

recovery_non_defol = read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/no_history/on_no_history_recovery_magnitude.csv")

# read in climate post
post_climate_no_history <-  read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/no_history/on_no_history_final_era5_clim.csv")

post_climate_history_1 <-  read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_history_final_era5_clim.csv")

post_climate_history_2 <-  read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_history_final_era5_clim_missing.csv")

# read in topography
topo <- read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_co-occurrences_topo.csv")

# merge recovery dataframes
recovery <- rbind(recovery_defol, recovery_non_defol) %>% 
  select(-c("X")) %>% 
  rename(recovery = Average.Recovery)

# merge post climate dataframes
post_climate_history <-  rbind(post_climate_history_1, post_climate_history_2)

# merge post climate dataframes
post_climate <-  rbind(post_climate_history, post_climate_no_history) 

# join recovery to history
history <- history %>% 
  left_join(recovery, by = "Fire_ID")

# join climate to history
history <- history %>% 
  left_join(post_climate, by = "Fire_ID")

# join fire centroids to history
history <- history %>% 
  left_join(centroids, by = "Fire_ID")

# join topography to history
history <- history %>% 
  left_join(topo, by = "Fire_ID")

# filter overlap percent to greater than 90
history_gt90 <- history %>% 
  filter(!(Max_Overlap_Percent <= 90 & history ==1 ))

#history greater than 90% as spatial object
h90.sf <- st_as_sf(history_gt90, coords = c("x", "y"), crs = 4326)

#sbw history
sbw <- read.csv("/home/goldma34/fire_insect_co-occurence/data/sbw-defol-data-v2.csv")

# clean history
history_gt90 <- history_gt90 %>% 
  left_join(sbw, by = "Fire_ID") %>% 
  select(-c(Time_Since_Defoliation, Cumulative_Years, defol)) %>% 
  rename(Time_Since_Defol = tsd) %>% 
  rename(Cumulative_Years_Defol = years_defol)


###############################################
#set windows of opp based on fleming et al 2002
##############################################


# set windows of opp
history_gt90 <- history_gt90 %>% 
  mutate(window_opp = case_when(Time_Since_Defol <= 2 & history ==1 ~ "1",
                                Time_Since_Defol >=3 & Time_Since_Defol<= 5 ~"2",
                                Time_Since_Defol >= 6 & Time_Since_Defol <=9 ~ "3",
                                Time_Since_Defol >= 10 ~ "4",
                                TRUE ~ "0"))


#Splitting the data for subclass
# keep non-defoliated options
hist_gt90_1 <- subset(history_gt90, window_opp == "0" | window_opp == "1")
hist_gt90_2  <- subset(history_gt90, window_opp == "0" | window_opp == "2")
hist_gt90_3 <- subset(history_gt90, window_opp == "0" | window_opp == "3")
hist_gt90_4 <- subset(history_gt90, window_opp == "0" | window_opp == "4")

#################################################
# set windows of opp with 0-2, 3-9, 10+
############################################

# set windows of opp
history_gt90 <- history_gt90 %>% 
  mutate(window_opp_2 = case_when(Time_Since_Defol <= 2 & history ==1 ~ "1",
                                Time_Since_Defol >=3 & Time_Since_Defol<= 9 ~"2",
                                Time_Since_Defol >= 10 ~ "3",
                                TRUE ~ "0"))


#Splitting the data for subclass
# keep non-defoliated options
hist_gt90_1_2 <- subset(history_gt90, window_opp_2 == "0" | window_opp_2 == "1")
hist_gt90_2_2  <- subset(history_gt90, window_opp_2 == "0" | window_opp_2 == "2")
hist_gt90_3_2 <- subset(history_gt90, window_opp_2 == "0" | window_opp_2 == "3")


#####################################
# matched to sf
#####################################


# load matched data
# matched data sf sev
m.data <- read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_sev_match_data.csv")
m.data.sev_sf <- st_as_sf(m.data, coords = c("x", "y"), crs = 4326)
# matched data sf rec
# Convert the data frame to an sf object
m.data.rec <- read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_rec_match_data.csv")
m.data.rec_sf <- st_as_sf(m.data.rec, coords = c("x", "y"), crs = 4326)