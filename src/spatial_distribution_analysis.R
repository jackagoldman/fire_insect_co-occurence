# Spatial Analysis

# load required packages
library(lme4)
library(MuMIn)
library(lmerTest)
library(nlme)
library(sp)
library(gstat)
library(ggplot2)
library(car)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(cowplot)


# Load your data 
source("/home/goldma34/fire_insect_co-occurence/src/load_data.R")  # Replace with actual script that loads hist_gt90_1

#load
source("/home/goldma34/fire_insect_co-occurence/src/find_best_model.R") 


# Load Ontario map data
ontario <- ne_states(country = "Canada", returnclass = "sf") %>%
  filter(name == "Ontario")


# Create the inset pie chart for m.data.sf
inset_data1 <- h90.sf %>%
  group_by(history) %>%
  summarise(count = n())

inset_pie1 <- ggplot(inset_data1, aes(x = "", y = count, fill = factor(history))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 3)  # Adjust the size here

# Create the inset pie chart for m.data.sf2
inset_data2 <- m.data.sf %>%
  group_by(history) %>%
  summarise(count = n())

inset_pie2 <- ggplot(inset_data2, aes(x = "", y = count, fill = factor(history))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar("y") +
  scale_fill_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                    name = NULL, 
                    labels = c("0" = "Non-Defoliated", "1" = "Defoliated")) +
  theme_void() +
  labs(x = NULL, y = NULL) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        legend.position = "none") +
  geom_text(aes(label = count), position = position_stack(vjust = 0.5), size = 3)  # Adjust the size here
# Create the first map with inset
map1 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = h90.sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                     name = NULL, 
                     labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                     guide = guide_legend(title.position = "top", title.hjust = 0.4)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = "right",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) +
  labs(title = "A)") +
  annotation_custom(ggplotGrob(inset_pie1), xmin = -78, xmax = -73, ymin = 42, ymax = 52)

# Create the second map with inset using m.data.sf2
map2 <- ggplot() +
  geom_sf(data = ontario, fill = "navajowhite1", color = "black") +
  geom_sf(data = m.data.sf, aes(color = factor(history)), size = 1) +
  scale_color_manual(values = c("0" = "#FF8C00A0", "1" = "#8B0000A0"), 
                     name = NULL, 
                     labels = c("0" = "Non-Defoliated", "1" = "Defoliated"),
                     guide = guide_legend(title.position = "top", title.hjust = 0.4)) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0, size = 16),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 10)
  ) +
  labs(title = "B)") +
  annotation_custom(ggplotGrob(inset_pie2), xmin = -78, xmax = -73, ymin = 42, ymax = 52)



# Extract the legend from one of the maps
legend <- cowplot::get_legend(map1)

# Combine the legend and the two maps
combined_plot <- plot_grid(
  legend,
  plot_grid(
    map1 + theme(legend.position = "none"),
    map2 + theme(legend.position = "none"),
    ncol = 2,
    rel_widths = c(1, 1)
  ),
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# Save the combined plot    
ggsave(
  filename = "/home/goldma34/fire_insect_co-occurence/data/figures/fig_1_map.png",
  plot = combined_plot,
  width = 12,
  height = 8,
  dpi = 300
)


# defol dataset
defol_only_sev <- subset(history_gt90, history == 1)
defol_only_rec <- subset(history_gt90, history == 1)

########################################
# Part 1: SEVERITY ======================
########################################
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
print(vif(model_gls))


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
