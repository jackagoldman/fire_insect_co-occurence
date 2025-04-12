# subgroup analysis

# load required packages
library(tidyverse)
library(DAPSm)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(cowplot)
library(MatchIt)
library(treatSens)
source("fire_insect_co-occurence/src/find_best_model.R")


# Required Data - prep -------------------------------------------------
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
  left_join(post_climate, by = 'Fire_ID')

# join fire centroids to history
history <- history %>% 
  left_join(centroids, by = "Fire_ID")

# join topography to history
history <- history %>% 
  left_join(topo, by = "Fire_ID")

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

# overview

# distribution of windows and observations
w1 <- hist_gt90_1 %>% group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",history ==0 ~"Non-Defoliated"))
w2 <- hist_gt90_2 %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated", history ==0 ~"Non-Defoliated"))
w3 <- hist_gt90_3 %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated", history ==0 ~"Non-Defoliated"))
w4 <- hist_gt90_4 %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",history ==0 ~"Non-Defoliated"))

w1 %>% filter(history == "Defoliated")
w2 %>% filter(history == "Defoliated")
w3 %>% filter(history == "Defoliated")
w4 %>% filter(history == "Defoliated")

# Part 1: Severity ==========================

# fit parameters
methods <- c("nearest", "optimal")
distances <- c("glm", "mahalanobis")
links <- c("logit", "probit")
m_orders <- c("random")
calipers <- c(0.1, 0.2, 0.3)
replace_list <- c(FALSE)
use_mahvars_list <- c(FALSE, TRUE)

## 1.1 baseline model =============
fit <- lm(rbr_w_offset  ~ history *(host_pct+ isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90),  data = m.data,weights = weights)

fit_att <- as.data.frame(marginaleffects::avg_comparisons(fit,
                                                          variables = "history",
                                                          vcov = ~subclass,
                                                          newdata = subset(history == 1))) %>% 
  mutate(Fires= "all")

## 1.2 subroup 1 ============
best_model_info <- find_best_model(hist_gt90_1, "severity",methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_1 <- best_model_info$best_model
all_results <- best_model_info$results

m.data.w1 <- match_data(best_model_1)


m.data.w1 %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",
                             history ==0 ~"Non-Defoliated"))

### 1.2.1 fit subgroup  ======
fit.w1 <- lm(rbr_w_offset  ~ history *(host_pct+ isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90),  data = m.data.w1,weights = weights)

fit.w1_att <- as.data.frame(marginaleffects::avg_comparisons(fit.w1,
                                                             variables = "history",
                                                             vcov = ~subclass,
                                                             newdata = subset(history == 1)))%>% 
  mutate(Fires= "0-2 years after defoliation")

### 1.2.3 sensitivity =========

vars.1<- dplyr::select(m.data.w1, c(host_pct, isi_90 , dc_90, dmc_90 , ffmc_90 , bui_90 , fwi_90 )) 
var_names<- names(vars.1)

X<- as.matrix(vars.1)
Y<- as.vector(m.data.w1$rbr_w_offset)
Z<- as.vector(m.data.w1$history)

# Define the sensitivity analysis model
sens_w1_sev <- treatSens(Y~ Z+X,
                         data = m.data.w1,
                         sensParam = "coef",  
                         trt.family = binomial(link = "probit"),  
                         grid.dim = c(10, 10),  
                         standardize = TRUE,
                         weights = "ATT", 
                         nsim = 100)

# Summarize the results
summary(sens_w1_sev)

sens_w1_sev$varnames<- c("Y","Z",var_names)

# save rData
saveRDS(sens_w1_sev, "fire_insect_co-occurence/data/results/sensitivity/severity_sensitivty_results_subgroup1.RDS")

sensPlot(sens_w1_sev, data.line=F, txtlab=T, which.txtlab=c(1,2,3,4,5,6,7))


## 1.3 subgroup 2 ================

best_model_info <- find_best_model(hist_gt90_2, "severity",methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_2 <- best_model_info$best_model
all_results <- best_model_info$results

m.data.w2 <- match_data(best_model_2)


m.data.w2 %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",
                             history ==0 ~"Non-Defoliated"))


### 1.3.1 fit subgroup  =======
fit.w2 <- lm(rbr_w_offset  ~ history *(host_pct+ isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90),  data = m.data.w2,weights = weights)

fit.w2_att <- as.data.frame(marginaleffects::avg_comparisons(fit.w2,
                                                             variables = "history",
                                                             vcov = ~subclass,
                                                             newdata = subset(history == 1)))%>% 
  mutate(Fires= "3-5 years after defoliation")

### 1.3.2 sensitivity ======

vars.2<- dplyr::select(m.data.w2, c(host_pct, isi_90 , dc_90, dmc_90 , ffmc_90 , bui_90 , fwi_90 )) 
var_names<- names(vars.2)

X<- as.matrix(vars.2)
Y<- as.vector(m.data.w2$rbr_w_offset)
Z<- as.vector(m.data.w2$history)

# Define the sensitivity analysis model
sens_w2_sev <- treatSens(Y~ Z+X,
                           data = m.data.w2,
                           sensParam = "coef",  
                           trt.family = binomial(link = "probit"),  
                           grid.dim = c(10, 10),  
                           standardize = TRUE,
                           weights = "ATT", 
                           nsim = 100)

# Summarize the results
summary(sens_w2_sev)

sens_w2_sev$varnames<- c("Y","Z",var_names)

# save rData
saveRDS(sens_w2_sev, "fire_insect_co-occurence/data/results/sensitivity/severity_sensitivty_results_subgroup2.RDS")

sensPlot(sens_w2_sev, data.line=F, txtlab=T, which.txtlab=c(1,2,3,4,5,6,7))

## subgroup 3 ==================
best_model_info_3 <- find_best_model(hist_gt90_3, "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_3 <- best_model_info_3$best_model
all_results_3 <- best_model_info_3$results



m.data.w3 <- match_data(best_model_3)


m.data.w3 %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",
                             history ==0 ~"Non-Defoliated"))

### 1.4.1 fit subgroup  =======
fit.w3 <- lm(rbr_w_offset  ~ history *(host_pct+ isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90),  data = m.data.w3,weights = weights)

fit.w3_att <- as.data.frame(marginaleffects::avg_comparisons(fit.w3,
                                                             variables = "history",
                                                             vcov = ~subclass,
                                                             newdata = subset(history == 1)))%>% 
  mutate(Fires= "6-9 years after defoliation")

### 1.4.2 sensitivity ======

vars.3<- dplyr::select(m.data.w3, c(host_pct, isi_90 , dc_90, dmc_90 , ffmc_90 , bui_90 , fwi_90 )) 
var_names<- names(vars.3)

X<- as.matrix(vars.3)
Y<- as.vector(m.data.w3$rbr_w_offset)
Z<- as.vector(m.data.w3$history)

# Define the sensitivity analysis model
sens_w3_sev <- treatSens(Y~ Z+X,
                         data = m.data.w3,
                         sensParam = "coef",  
                         trt.family = binomial(link = "probit"),  
                         grid.dim = c(10, 10),  
                         standardize = TRUE,
                         weights = "ATT", 
                         nsim = 100)

# Summarize the results
summary(sens_w3_sev)

sens_w3_sev$varnames<- c("Y","Z",var_names)

# save rData
saveRDS(sens_w3_sev, "fire_insect_co-occurence/data/results/sensitivity/severity_sensitivty_results_subgroup3.RDS")
sensPlot(sens_w3_sev, data.line=F, txtlab=T, which.txtlab=c(1,2,3,4,5,6,7))


## subgroup 4 ==================
best_model_info_4 <- find_best_model(hist_gt90_4, "severity", methods, distances, links, m_orders, calipers, replace_list, use_mahvars_list)
best_model_4 <- best_model_info$best_model
all_results_4 <- best_model_info$results



m.data.w4 <- match_data(best_model_4)


m.data.w4 %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",
                             history ==0 ~"Non-Defoliated"))

### 1.4.1 fit subgroup  =======
fit.w4 <- lm(rbr_w_offset  ~ history *(host_pct+ isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90),  data = m.data.w4,weights = weights)

fit.w4_att <- as.data.frame(marginaleffects::avg_comparisons(fit.w4,
                                                             variables = "history",
                                                             vcov = ~subclass,
                                                             newdata = subset(history == 1)))%>% 
  mutate(Fires= "10-15 years after defoliation")

### 1.5.2 sensitivity ======

vars.4<- dplyr::select(m.data.w4, c(host_pct, isi_90 , dc_90, dmc_90 , ffmc_90 , bui_90 , fwi_90 )) 
var_names<- names(vars.4)

X<- as.matrix(vars.4)
Y<- as.vector(m.data.w4$rbr_w_offset)
Z<- as.vector(m.data.w4$history)

# Define the sensitivity analysis model
sens_w4_sev <- treatSens(Y~ Z+X,
                         data = m.data.w4,
                         sensParam = "coef",  
                         trt.family = binomial(link = "probit"),  
                         grid.dim = c(10, 10),  
                         standardize = TRUE,
                         weights = "ATT", 
                         nsim = 100)

# Summarize the results
summary(sens_w4_sev)

sens_w4_sev$varnames<- c("Y","Z",var_names)

# save rData
saveRDS(sens_w4_sev, "fire_insect_co-occurence/data/results/sensitivity/severity_sensitivty_results_subgroup4.RDS")
sensPlot(sens_w4_sev, data.line=F, txtlab=T, which.txtlab=c(1,2,3,4,5,6,7))
