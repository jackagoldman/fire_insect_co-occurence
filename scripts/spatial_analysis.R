# Spatial Analysis

# load required packages
library(lme4)
library(MuMIn)
library(lmerTest)
library(nlme)
library(sp)
library(gstat)
library(ggplot2)
source("fire_insect_co-occurence/src/find_best_model.R")


# Prep data ===========
# read in history data
history = read.csv("/home/goldma34/fire_insect_co-occurence/data/outputs/on/on_defoliation_history_wx.csv")

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
  left_join(post_climate, by = 'Fire_ID')

# join fire centroids to history
history <- history %>% 
  left_join(centroids, by = "Fire_ID")

#
history_gt90 <- history %>% 
  filter(!(Max_Overlap_Percent <= 90 & history ==1 ))

#
h90.sf <- st_as_sf(history_gt90, coords = c("x", "y"), crs = 4326)

#sbw history
sbw <- read.csv("/home/goldma34/fire_insect_co-occurence/data/sbw-defol-data-v2.csv")

# clean history
history_gt90 <- history_gt90 %>% 
  left_join(sbw, by = "Fire_ID") %>% 
  select(-c(Time_Since_Defoliation, Cumulative_Years, defol)) %>% 
  rename(Time_Since_Defol = tsd) %>% 
  rename(Cumulative_Years_Defol = years_defol)


# set windows of opp
history_gt90 <- history_gt90 %>% 
  mutate(window_opp = case_when(Time_Since_Defol <= 3 & history ==1 ~ "1",
                                Time_Since_Defol >=4 & Time_Since_Defol<= 6 ~"2",
                                Time_Since_Defol >= 7 & Time_Since_Defol <=9 ~ "3",
                                Time_Since_Defol >= 10 ~ "4",
                                TRUE ~ "0"))


#Splitting the data for subclass
# keep non-defoliated options
hist_gt90_1 <- subset(history_gt90, window_opp == "0" | window_opp == "1")
hist_gt90_2  <- subset(history_gt90, window_opp == "0" | window_opp == "2")
hist_gt90_3 <- subset(history_gt90, window_opp == "0" | window_opp == "3")
hist_gt90_4 <- subset(history_gt90, window_opp == "0" | window_opp == "4")


# defol dataset
defol_only_sev <- subset(history_gt90, history == 1)
defol_only_rec <- subset(history_gt90, history == 1)


# Part 1: SEVERITY ======================
# Define the model formula
formula <- rbr_w_offset ~host_pct +  Cumulative_Years_Defol:window_opp + isi_90 + 
  dc_90 + dmc_90 + ffmc_90 + bui_90+ fwi_90 + mean_tri+  x + y 

# Fit the model using nlme with spatial correlation, where correlation decreases with distance
# gls
# Fit the model using nlme with spatial correlation, where correlation decreases with distance
model_gls <- gls(formula,
                 correlation = corGaus(form = ~ x + y), 
                 data = defol_only_sev,
                 control = list(singular.ok = TRUE))

# Check for multicollinearity
library(car)
print(vif(model_gls))

model_gls <- gls(formula,
                 correlation = corGaus(form = ~ x + y), 
                 data = defol_only)


test_data <- droplevels(defol_only)

# Residuals vs. Fitted Values Plot
resid_plot_severity <- plot(fitted(model_gls), resid(model_gls)) +abline(h = 0, col = "red") # fit is good

# Display the model summary
summary(model_gls)

# r2
rsquared.gls(model_gls, formula)

# extract residuals
defol_only_sev$residuals_gls <- residuals(model_gls)

colnames(defol_only_sev) <- make.names(colnames(defol_only_sev))

# get coords
coordinates(defol_only_sev) <- ~ x + y

# Compute the variogram
variogram_gls <- variogram(residuals_gls ~ 1, data = defol_only_sev)

# Plot the variogram
var_plot <- plot(variogram_gls, main = "Variogram of Residuals (Gaussian)")

# confirm defol_only as dataframe
df <- as.data.frame(defol_only_sev)

# Create a spatial plot of residuals
ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df, aes(x = x, y = y, color = residuals_gls)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()

# Part 2: RECOVERY ==================

# Define the model formula
formula <- recovery ~host_pct +  Cumulative_Years_Defol:window_opp + rbr_w_offset +
  mean_temperature + sum_precipitation_mm + mean_tri+ x + y 

# Fit the model using nlme with spatial correlation, where correlation decreases with distance
# gls

model_gls_rec <- gls(formula,
                     correlation = corGaus(form = ~ x + y), 
                     data = defol_only_rec)

# Residuals vs. Fitted Values Plot
residual_plot_rec <- plot(fitted(model_gls_rec), resid(model_gls_rec)) +abline(h = 0, col = "red") # fit is good

#model_glsI <- update(model_gls, weights = varIdent(form = ~ 1 | Time_Since_Defol))


# Display the model summary
summary(model_gls_rec)

# r2
rsquared.gls(model_gls_rec, formula)

# extract residuals
defol_only_rec$residuals_gls_rec <- residuals(model_gls_rec)

colnames(defol_only_rec) <- make.names(colnames(defol_only_rec))

# get coords
coordinates(defol_only_rec) <- ~ x + y

# Compute the variogram
variogram_gls_rec <- variogram(residuals_gls_rec ~ 1, data = defol_only_rec)

# Plot the variogram
var_plot_rec <- plot(variogram_gls_rec, main = "Variogram of Residuals (Gaussian)")

# confirm defol_only as dataframe
df_rec <- as.data.frame(defol_only_rec)

# Create a spatial plot of residuals
ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_point(data = df_rec, aes(x = x, y = y, color = residuals_gls_rec)) +
  scale_color_gradient2(low = "#88CCEEA0", mid = "white", high = "#8B0000A0", midpoint = 0) +
  labs(title = NULL, x = "Longitude", y = "Latitude", color = "Residuals") +
  theme_bw()
