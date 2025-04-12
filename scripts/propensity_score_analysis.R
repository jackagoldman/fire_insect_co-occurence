# Propensity score analysis 


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
source("src/find_best_model.R")


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





# Part 1: SEVERITY ==============================

## Propensity score matching for full severity Data ---------------------------

# Full matching on a probit PS
m.out2 <- matchit(history ~  host_pct +isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri ,
                  data = history_gt90, 
                  method = "nearest",
                  distance = "glm",
                  link = "probit",
                  m.order = "random",
                  ### set a maximum distance for the matches
                  caliper = 0.25)



m.data <- match_data(m.out2)


m.data %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",
                             history ==0 ~"Non-Defoliated"))

# Convert the data frame to an sf object
m.data.sf <- st_as_sf(m.data, coords = c("x", "y"), crs = 4326)


## Assess balance for severity ======================
v <- data.frame(old = c("host_pct", "isi_90", "dc_90", 
"dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri"),
new = c("Host Species Percentage", "Intial Spread Index", "Drought Code", 
        "Duff Moisture Code", "Fine Fuel Moisture Code", "Build Up Index", "Fire Weather Index", "Mean Terrain Ruggedness Index"))

cobalt::love.plot(m.out2, stats = c("mean.diffs"), 
                  threshold = c(m = .25), 
                  binary = "std",
                  abs = TRUE,
                  var.order = "unadjusted",
                  var.names = v,
                  limits = c(0, 1),
                  grid = FALSE,
                  wrap = 20,
                  sample.names = c("Unmatched", "Matched"),
                  position = "top",
                  shapes = c("circle", "triangle"),
                  colors = c("#FF8C00A0", "#8B0000A0"))


host_plot <- cobalt::bal.plot(m.out2, 
                              var.name = "host_pct", 
                              which = "both",
                              colors = c("#FF8C00A0", "#8B0000A0"),
                              sample.names = c("Unmatched", "Matched")) +scale_fill_manual(
                                values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                                labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Host Species Percentage") + ggtitle(NULL) +theme_bw() + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

isi_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "isi_90", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Initial Spread Index") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

dc_plot <- cobalt::bal.plot(m.out2, 
                            var.name = "dc_90", 
                            which = "both", 
                            colors = c("#FF8C00A0", "#8B0000A0"),
                            sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                              values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                              labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Drought Code") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

dmc_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "dmc_90", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Duff Moisture Code") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

ffmc_plot <- cobalt::bal.plot(m.out2, 
                              var.name = "ffmc_90", 
                              which = "both", 
                              colors = c("#FF8C00A0", "#8B0000A0"),
                              sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                                values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                                labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Fine Fuel Moisture Code") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

bui_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "bui_90", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Build Up Index") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

fwi_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "fwi_90", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Fire Weather Index") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

tri_plot <- cobalt::bal.plot(m.out2, 
                             var.name = "mean_tri", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Fire Weather Index") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

# Arrange the plots in a 2-column layout
plot_grid(host_plot, isi_plot, dc_plot, dmc_plot, ffmc_plot, bui_plot,fwi_plot, ncol = 1, align = 'v', rel_heights = c(1, 1, 1, 1, 1, 1))


## Severity average treatment effect estimate =====================

# fit model
fit_sev <- lm(rbr_w_offset  ~ history + host_pct+ isi_90 + dc_90+ dmc_90 + ffmc_90 + bui_90 + fwi_90 + mean_tri,  data = m.data,weights = weights)

# avg_comparisons for treatment effect
fit_att_sev <- as.data.frame(marginaleffects::avg_comparisons(fit_sev,
                                                          variables = "history",
                                                          vcov = ~subclass,
                                                          newdata = subset(history == 1))) %>%  mutate(Fires= "all")

# results
fit_att_sev

#rsquared
r.squaredGLMM(fit_sev)

## Sensitivity analysis for severity #####

# benchmark covariates
covariates <- c("host_pct", "isi_90", "dc_90", "dmc_90", "ffmc_90", "bui_90", "fwi_90", "mean_tri")

# run sensitivity analysis
sev.sensitivity <- sensemakr(model = fit_sev, 
                                treatment = "history",
                                benchmark_covariates = covariates,
                                alpha = 0.05, 
                                reduce = TRUE)

#outcome
summary(sev.sensitivity)

#The results suggest that the 'history' variable has a moderate coefficient estimate with a relatively 
#low standard error, leading to a high t-value and indicating that the effect is statistically significant. 
#The sensitivity statistics imply that the observed effect is reasonably robust and would 
#require relatively strong unobserved confounders to nullify or make the effect statistically insignificant.

#plot bias contour of point estimate
plot(sev.sensitivity)

# plot extreme scenario
plot(sev.sensitivity, type = "extreme")


# Part 2: RECOVERY ==============================

## Propensity score matching for full recovery Data --------------------------

# Full matching on a probit PS
m.out.rec <- matchit(history ~  host_pct + rbr_w_offset + mean_temperature + sum_precipitation_mm + mean_tri,
                     data = history_gt90,
                     method = "nearest",
                     distance = "glm",
                     link = "probit",
                     m.order = "random",
                     ### set a maximum distance for the matches
                     caliper = 0.25)



m.data.rec <- match_data(m.out.rec)


m.data.rec %>% 
  group_by(history) %>% 
  summarise("Number of Fires" = n()) %>% 
  mutate(history = case_when(history ==1 ~ "Defoliated",
                             history ==0 ~"Non-Defoliated"))

# Convert the data frame to an sf object
m.data.rec.sf <- st_as_sf(m.data.rec, coords = c("x", "y"), crs = 4326)


## Assess balance for recovery ==============
v.r <- data.frame(old = c("host_pct", "rbr_w_offset", "mean_temperature", "sum_precipitation_mm"),
                new = c("Host Species Percentage", "burn severity", "post-fire mean temperature", 
                        "post-fire total precipitation (mm)"))

cobalt::love.plot(m.out.rec, stats = c("mean.diffs"), 
                  threshold = c(m = .25), 
                  binary = "std",
                  abs = TRUE,
                  var.order = "unadjusted",
                  var.names = v.r,
                  limits = c(0, 1),
                  grid = FALSE,
                  wrap = 20,
                  sample.names = c("Unmatched", "Matched"),
                  position = "top",
                  shapes = c("circle", "triangle"),
                  colors = c("#FF8C00A0", "#8B0000A0"))


host_plot <- cobalt::bal.plot(m.out.rec, 
                              var.name = "host_pct", 
                              which = "both",
                              colors = c("#FF8C00A0", "#8B0000A0"),
                              sample.names = c("Unmatched", "Matched")) +scale_fill_manual(
                                values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                                labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Host Species Percentage") + ggtitle(NULL) +theme_bw() + theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

rbr_plot <- cobalt::bal.plot(m.out.rec, 
                             var.name = "rbr_w_offset", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Burn Severity (RBR)") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

temp_plot <- cobalt::bal.plot(m.out.rec, 
                            var.name = "mean_temperature", 
                            which = "both", 
                            colors = c("#FF8C00A0", "#8B0000A0"),
                            sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                              values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                              labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Mean Temperature (°C)") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

ppt_plot <- cobalt::bal.plot(m.out.rec, 
                             var.name = "sum_precipitation_mm", 
                             which = "both", 
                             colors = c("#FF8C00A0", "#8B0000A0"),
                             sample.names = c("Unmatched", "Matched"))+ scale_fill_manual(
                               values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"),
                               labels = c("0" = "Non-Defoliated", "1" = "Defoliated"))+xlab("Total Precipitation (mm)") + ggtitle(NULL) +theme_bw()+ theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))


# Arrange the plots in a 2-column layout
plot_grid(host_plot, rbr_plot, temp_plot, ppt_plot, ncol = 1, align = 'v', rel_heights = c(1, 1, 1, 1, 1, 1))


## Recovery average treatment estimate ========
m.data.rec_test  <- m.data.rec %>% 
  mutate(history = case_when(history == 1 ~ "Defoliated",
                             history == 0 ~ "Non-Defoliated"))

fit_rec_baseline <- lm(recovery  ~ history + host_pct+rbr_w_offset + mean_temperature + sum_precipitation_mm + mean_tri,  data = m.data.rec,weights = weights)

fit_att_rec <- as.data.frame(marginaleffects::avg_comparisons(fit_rec_baseline,
                                                              variables = "history",
                                                              vcov = ~subclass,
                                                              newdata = subset(history == 1))) %>% 
  mutate(Fires= "all")

# summary
fit_att_rec

#r squared
r.squaredGLMM(fit_rec_baseline)

# effects plot
p <- marginaleffects::plot_predictions(fit_rec_baseline, condition = c("history"), draw = FALSE) 

p <- p %>% 
  mutate(Fires = "All")

plot1 <- ggplot(p, aes(x = history, y = estimate, color = history)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) +
  labs(y = "Recovery Magnitude (%)", x = "Defoliation History", color = "History") +
  scale_x_discrete(labels = c("Non-Defoliated", "Defoliated")) +
  scale_color_manual(values = c("Defoliated" = "#FF8C00A0", "Non-Defoliated" = "#8B0000A0")) +
  theme_bw()+
  theme(legend.position = "none")

##  sensitivity analysis for recovery  -----------

# benchmark covariates
covariates_rec <- c("host_pct", "rbr_w_offset", "mean_temperature", "sum_precipitation_mm", "mean_tri")

# run sensitivity analysis
rec.sensitivity <- sensemakr(model = fit_rec_baseline, 
                             treatment = "history",
                             benchmark_covariates = covariates_rec,
                             alpha = 0.05, 
                             reduce = TRUE)

#outcome
summary(rec.sensitivity)

#The results suggest that the 'history' variable has a moderate coefficient estimate with a relatively 
#low standard error, leading to a high t-value and indicating that the effect is statistically significant. 
#The sensitivity statistics imply that the observed effect is reasonably robust and would require relatively 
#strong unobserved confounders to nullify or make the effect statistically insignificant.
#Overall, these results indicate that 'history' has a significant and robust effect on the outcome.

#plot bias contour of point estimate
plot(rec.sensitivity)

# plot extreme scenario
plot(rec.sensitivity, type = "extreme")


