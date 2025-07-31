# Import Libraries
library(tidyverse)
library(readxl)
library(sf)
library(tmap)
library(RColorBrewer)
library(car)
library(units)
library(lmtest)
library(sandwich)
library(ggplot2)
library(terra)
library(exactextractr)
library(installr)
library(e1071)
library(elevatr)
library(GWmodel)

# Set the working directory
setwd("C:/Users/andre/OneDrive - Georgia Institute of Technology/Georgia Tech/Classes/Summer 2025/CP 6596/")
getwd()

# Set Parameters
stateOfInterest <- c("GA")
allStates <- c("US")
projected_crs <- 26916  # NAD83 / UTM zone 16N

# Load in Data
# County Polygon shapefile
ga_counties = st_read("Data/tl_2020_13_county20.shp")
# City Polygon shapefile
ga_cities = st_read("Data/Cities_2019_TIGER/Cities_2019_TIGER.shp")
# Watersheds
watershed_boundaries = st_read("Data/WBD_03_HU2_Shape/Shape/WBDLine.shp")
watersheds = st_read("Data/WBD_03_HU2_Shape/Shape/WBDHU12.shp")

# Process counties and cities to create municipalities map
#st_crs(ga_cities)
#st_crs(ga_counties)
# Update CRS
cities <- st_transform(ga_cities, st_crs(ga_counties))
#st_crs(cities)
# Dissolve all cities into one geometry for clipping purposes
incorporated_union <- st_union(cities)

# Clip counties to remove incorporated areas (i.e., get unincorporated areas)
unincorporated <- st_difference(ga_counties, incorporated_union)

# Add a new attribute to each dataset
cities <- cities %>% mutate(status = "Incorporated")
unincorporated <- unincorporated %>% mutate(status = "Unincorporated")

# Combine them into one shapefile
municipalities <- rbind(
  select(cities, status, geometry),
  select(unincorporated, status, geometry)
)

watershed_boundaries <- st_transform(watershed_boundaries, projected_crs)
municipalities <- st_transform(municipalities, projected_crs)
watersheds <- st_transform(watersheds, projected_crs)

watershed_boundaries <- st_filter(watershed_boundaries, municipalities, .predicate = st_intersects)
watersheds <- st_filter(watersheds, municipalities, .predicate = st_intersects)

# Write to a new shapefile
# st_write(municipalities, "incorporated_unincorporated_combined.shp")

# Generate boundaries from municipalities
boundary_lines <- st_cast(st_boundary(municipalities), "LINESTRING")
# Remove invalid
boundary_lines <- boundary_lines[st_is_valid(boundary_lines), ]
# Remove lines with length < 1m
boundary_lines <- boundary_lines[as.numeric(st_length(boundary_lines)) > 1, ]

# Precompute municipal areas
municipalities$area_m2 <- as.numeric(st_area(municipalities))


# Prepare results dataframe
results <- data.frame(
  watershed = character(),
  intersection_count = numeric(),
  total_boundary_length_km = numeric(),
  weighted_boundary_length_km = numeric(),
  dominance_ratio = numeric(),
  area_km2 = numeric(),
  id = character(),
  stringsAsFactors = FALSE
)

# Loop through watersheds
for (i in 1:nrow(watersheds)) {
  ws <- watersheds[i, ]
  ws_name <- ws$name
  ws_geom <- ws$geometry
  
  # Check if any intersection
  if (any(st_intersects(ws, municipalities, sparse = FALSE))) {
    
    # Intersect once
    lines_in_ws <- st_intersection(boundary_lines, ws)
    
    lines_in_ws <- st_transform(lines_in_ws, projected_crs)
    
    # Cast to LINESTRING
    lines_in_ws <- st_cast(lines_in_ws, "LINESTRING")
    line_lengths <- st_length(lines_in_ws)
    
    if (nrow(lines_in_ws) > 0) {
      
      # Total boundary length (dissolve overlapping lines)
      lines_union <- st_union(lines_in_ws)
      total_length_km <- sum(st_length(lines_union)) / 1000
      lines_union <- st_cast(lines_union, "LINESTRING")
      print("Total length (km):")
      print(total_length_km)
      
      # Weighted length (based on distance to watershed centroid)
      ws_centroid <- st_centroid(ws)
      
      # Sample points along lines
      #points <- st_line_sample(lines_union, density = 1) # one point per 1m
      #points_sf <- st_sf(geometry = points)
      #points_list <- lapply(1:length(lines_in_ws), function(i) {
      #  n_points <- as.integer(line_lengths[i])  # 1 point per meter
      #  print(n_points)
      #  if (!is.na(n_points) && n_points > 0) {
      #    st_line_sample(lines_in_ws[i, ], n = n_points, type = "regular")
      #  } else {
      #    st_sfc()
      #  }
      #})
      #points_sf <- do.call(c, points)
      
      
      # Densify lines — segmentize so every segment is ~1m (adjust as needed)
      densified_lines <- st_segmentize(lines_in_ws, dfMaxLength = 1)
      
      print("Densified lines:")
      print(length(densified_lines))
      
      # Extract points along the densified lines
      points <- st_cast(densified_lines, "POINT")
      points <- st_transform(points, projected_crs)
      
      
      print("Number of points:")
      print(length(points_sf))
      
      # Compute max distance from centroid to watershed boundary
      ws_boundary_pts <- st_cast(st_boundary(ws), "LINESTRING") %>% 
        st_segmentize(dfMaxLength = 10) %>% 
        st_cast("POINT")      
      max_dist <- max(st_distance(ws_centroid, ws_boundary_points))
      # Compute distances to centroid and weights
      distances <- st_distance(points, ws_centroid)
      
      weights <- 0.5 + (as.numeric(distances) / as.numeric(max_dist))  # closer points get higher weight
      
      print("Distances:")
      print(head(distances))
      
      print("Accumulated Weighted points:")
      print(summary(weights))
      
      # Approximate weighted length: total weights in km
      if (length(weights) > 0) {
        weighted_length_km <- (sum(weights) / 1000)
      } else {
        weighted_length_km <- 0
      }
    } else {
      total_length_km <- 0
      weighted_length_km <- 0
    }
    
    # Dominance ratio: find area per muni in this watershed
    munis_in_ws <- st_intersection(municipalities, ws)
    munis_in_ws$area <- as.numeric(st_area(munis_in_ws))
    
    #print(munis_in_ws)
    
    # Intersection count: unique muni names intersecting
    count <- nrow(munis_in_ws)
    
    #print(nrow(munis_in_ws))
    if (nrow(munis_in_ws) > 1) {
      muni_areas <- sort(unlist(munis_in_ws$area_m2),decreasing=TRUE)
      #print(nrow(muni_areas))
      #print(muni_areas)
      top_area <- unlist(muni_areas[1])
      #print("Top area:")
      #print(top_area)
      smaller_muni_areas <- muni_areas[-1]
      other_areas <- sum(unlist(smaller_muni_areas))
      #print("Other areas:")
      #print(other_areas)
      dominance_ratio <- (top_area - other_areas)/(top_area + other_areas)
    } else {
      dominance_ratio <- 1
    }
    
  } else {
    total_length_km <- 0
    weighted_length_km <- 0
    count <- 0
    dominance_ratio <- 0
  }
  
  # Store result
  results <- rbind(results, data.frame(
    watershed = ws_name,
    intersection_count = count,
    total_boundary_length_km = total_length_km,
    weighted_boundary_length_km = weighted_length_km,
    dominance_ratio = dominance_ratio,
    area_km2 = ws$areasqkm,
    id = ws$tnmid,
    stringsAsFactors = FALSE
  ))
  
  print(paste0(ws_name, ": ", count, " munis, ", round(total_length_km, 2), " km total, ", 
               round(weighted_length_km, 2), " km weighted, dominance ratio = ", round(dominance_ratio, 2)))
}

# Normalize by watershed area
results_normalized <- results %>%
  mutate(tnmid = id) %>%
  mutate(normalized_intersections = intersection_count / area_km2) %>%
  mutate(normalized_total_length = total_boundary_length_km / area_km2) %>%
  mutate(normalized_weighted_length = weighted_boundary_length_km / area_km2)


watershed_statistics <- results_normalized

# Environmental Justice Data

ejscreen = st_read("Data/EJScreen.shp")
view(ejscreen)
ej <- ejscreen %>% filter(ST_ABBREV %in% stateOfInterest)# %>% select()
#drinking_water <- ej %>% select(ID, DWATER, D2_DWATER, D5_DWATER)
#view(drinking_water)

#print(names(ej))
#str(ej$DWATER)

# Check and align CRS
if (st_crs(ej) != st_crs(watersheds)) {
  ej <- st_transform(ej, st_crs(results))
}

# Spatially join: intersect block groups with watersheds
# This will split block groups where they cross watersheds
ej_watersheds <- st_intersection(ej, watersheds)

# Calculate area of intersection for weighting (if you don't have population)
ej_watersheds$area_overlap <- st_area(ej_watersheds)

# Aggregate: compute area-weighted average DWATER per watershed
watershed_dwater <- ej_watersheds %>%
  group_by(tnmid) %>%
  summarise(
    DWATER_mean = sum(as.numeric(DWATER) * as.numeric(area_overlap), na.rm = TRUE) / sum(as.numeric(area_overlap), na.rm = TRUE),
    pop_mean = sum(as.numeric(ACSTOTPOP) * as.numeric(area_overlap), na.rm = TRUE) / sum(as.numeric(area_overlap), na.rm = TRUE),
    wastewater_discharge_mean = sum(as.numeric(PWDIS) * as.numeric(area_overlap), na.rm = TRUE) / sum(as.numeric(area_overlap), na.rm = TRUE),
    lowincome_mean = sum(as.numeric(LOWINCPCT) * as.numeric(area_overlap), na.rm = TRUE) / sum(as.numeric(area_overlap), na.rm = TRUE),
    .groups = "drop"
  )

ws_dwater_df = st_drop_geometry(watershed_dwater)

# Join aggregated result back to original watersheds layer
watershed_statistics <- left_join(watershed_statistics, ws_dwater_df, by = "tnmid")

wildfires <- st_read("Data/Wildfires_16N.shp")

st_crs(watersheds)
st_crs(wildfires)
if (st_crs(wildfires) != st_crs(watersheds)) {
  wildfires <- st_transform(wildfires, st_crs(watersheds))
}

unique(wildfires$DROUGHT)

# Filter drought-caused wildfires
drought_fires <- wildfires %>%
  filter(DROUGHT >= 100)

unique(drought_fires$SIZE)

# Assign numeric severity score based on size class
size_weights <- c("A" = 0.1, "B" = 5, "C" = 50, "D" = 150, "E" = 500, "F" = 2500, "G" = 5000)

drought_fires <- drought_fires %>%
  mutate(
    severity_score = size_weights[SIZE]
  )

fires <- wildfires %>%
  mutate(
    severity_score = size_weights[SIZE]
  )

# Spatial join with watersheds
fires_with_ws <- st_join(fires, watersheds, join = st_intersects, left = FALSE)

# Aggregate by watershed
fire_summary <- fires_with_ws %>%
  st_drop_geometry() %>%
  group_by(tnmid) %>%
  summarise(
    #drought_fire_count = n(),
    fire_count = n(),
    total_severity = sum(severity_score, na.rm = TRUE),
    mean_severity = mean(severity_score, na.rm = TRUE),
    max_severity = max(severity_score, na.rm = TRUE)
  )

# Join back to watersheds
watershed_statistics <- left_join(watershed_statistics, fire_summary, by = "tnmid")

# Compute per km²
watershed_statistics <- watershed_statistics %>%
  mutate(
    severity_per_km2 = total_severity / area_km2,
    fires_per_km2 = fire_count / area_km2
  )

# Write results to CSV# Write results to CSVwatershed
write.csv(watershed_statistics, "watershed_stats_week10.csv", row.names = FALSE)

droughts <- read.csv("Data/ncei_noaa_pdsi_ga.csv") %>% mutate(Name = str_remove(Name, " County$"))

drought_counties <- ga_counties  %>% 
  left_join(droughts, by = c("NAME20" = "Name"))

st_crs(watersheds)
st_crs(drought_counties)
if (st_crs(drought_counties) != st_crs(watersheds)) {
  drought_counties <- st_transform(drought_counties, st_crs(watersheds))
}

summary(drought_counties$X1901.2000.Mean)

drought_watersheds <- st_intersection(drought_counties, watersheds)

# Calculate area of intersection for weighting (if you don't have population)
drought_watersheds$area_overlap <- st_area(drought_watersheds)

# Aggregate: compute area-weighted average drought per watershed
watershed_droughts <- drought_watersheds %>%
  group_by(tnmid) %>%   
  summarise(
    DROUGHT_mean = sum(as.numeric(X1901.2000.Mean) * as.numeric(area_overlap), na.rm = TRUE) / sum(as.numeric(area_overlap), na.rm = TRUE),
    AWATER_mean = sum(as.numeric(AWATER20) * as.numeric(area_overlap), na.rm = TRUE) / sum(as.numeric(area_overlap), na.rm = TRUE),
    .groups = "drop"
  )

ws_drought_df = st_drop_geometry(watershed_droughts)

# Join aggregated result back to original watersheds layer
watershed_statistics <- left_join(watershed_statistics, ws_drought_df, by = "tnmid")

write.csv(watershed_statistics_drought, "watershed_stats_drought.csv", row.names = FALSE)

floods <- st_read("Data/NRI_Shapefile_CensusTracts/NRI_Shapefile_CensusTracts.shp")

st_crs(watersheds)
st_crs(floods)
if (st_crs(floods) != st_crs(watersheds)) {
  floods <- st_transform(floods, st_crs(watersheds))
}

ga_floods <- floods %>% filter(STATE == "Georgia")

floods_watersheds <- st_intersection(ga_floods, watersheds)

# Calculate area of intersection for weighting (if you don't have population)
floods_watersheds$area_overlap <- st_area(floods_watersheds)

# Aggregate: compute area-weighted average floods per watershed
watershed_floods <- floods_watersheds %>%
  group_by(tnmid) %>%   
  summarise(
    FLOOD_mean = sum(as.numeric(RISK_SCORE) * as.numeric(area_overlap), na.rm = TRUE) / sum(as.numeric(area_overlap), na.rm = TRUE),
    .groups = "drop"
  )

ws_floods_df = st_drop_geometry(watershed_floods)

# Join aggregated result back to original watersheds layer
watershed_statistics <- left_join(watershed_statistics, ws_floods_df, by = "tnmid")

write.csv(watershed_statistics, "watershed_stats_floods.csv", row.names = FALSE)

if (1) {
  #packageVersion("tmap")
  tmap_mode("view")
  
  # Create interactive thematic map
  Figure <- tm_shape(ej) +
    tm_fill(col = "DWATER",
            fill.scale = tm_scale_continuous(values = "brewer.reds")) +
    tm_layout(legend.outside = TRUE) +
    tm_title("Drinking Water")
  
  # View TMap
  Figure
  # tm_view(Figure)
}

# ADD ELEVATION DATA
# Get a gridded DEM clipped to watershed extent
httr::set_config(httr::timeout(6000))  # 60 seconds
dem10m <- get_elev_raster(watersheds, z = 10)  # z 10 >≈ 10‑30 m resolution
dem10m

# Extract mean elevation
ws_elev <- terra::extract(rast(dem10m), vect(watersheds), fun = mean, na.rm = TRUE)
watershed_statistics$mean_elevation <- ws_elev[,2]

write.csv(watershed_statistics, "watershed_stats_elevation.csv", row.names = FALSE)

# Consolidate independent variables
# watershed_statistics <- results_normalized %>% mutate(name = watershed)

watershed_statistics <- watershed_statistics %>% mutate(boundary_centrality_score = (weighted_boundary_length_km - total_boundary_length_km) / weighted_boundary_length_km)
watershed_statistics <- watershed_statistics %>% mutate(boundary_centrality_score = replace_na(boundary_centrality_score, -1))
summary(watershed_statistics$boundary_centrality_score)

watershed_statistics <- watershed_statistics %>% mutate(boundary_centrality_index = (weighted_boundary_length_km - total_boundary_length_km) / area_km2)
summary(watershed_statistics$boundary_centrality_index)

watershed_statistics <- watershed_statistics %>% mutate(boundary_centrality_ratio = weighted_boundary_length_km / total_boundary_length_km)
watershed_statistics <- watershed_statistics %>% mutate(boundary_centrality_ratio = replace_na(boundary_centrality_ratio, 0))
summary(watershed_statistics$boundary_centrality_ratio)

watershed_statistics <- watershed_statistics %>% mutate(boundary_centrality = weighted_boundary_length_km - total_boundary_length_km)#/ area)
summary(watershed_statistics$boundary_centrality)

watershed_statistics <- watershed_statistics %>% mutate(combined_boundary_length = weighted_boundary_length_km + total_boundary_length_km)#/ area)
summary(watershed_statistics$combined_boundary_length)

watershed_statistics <- watershed_statistics %>% mutate(normalized_combined_boundary_length = weighted_boundary_length_km + total_boundary_length_km/ area_km2)
summary(watershed_statistics$normalized_combined_boundary_length)

summary(watershed_statistics$intersection_count)
watershed_statistics <- watershed_statistics %>% mutate(watershed_dividedness = boundary_centrality_index + intersection_count/6 + dominance_ratio/2)
summary(watershed_statistics$watershed_dividedness)

summary(watershed_statistics$DWATER_mean)
summary(watershed_statistics$DROUGHT_mean)
summary(watershed_statistics$FLOOD_mean)
summary(watershed_statistics$wastewater_discharge_mean)
watershed_statistics <- watershed_statistics %>% mutate(watershed_burden = DWATER_mean + wastewater_discharge_mean/100 + FLOOD_mean + abs(DROUGHT_mean * 100))
summary(watershed_statistics$watershed_burden)
summary(log(watershed_statistics$watershed_burden))


# Find watersheds with NA flood mean values
watersheds_to_remove <- watershed_statistics %>%
  filter(is.na(FLOOD_mean)) %>%
  pull(watershed)  

#watershed_statistics <- watershed_statistics %>% select(-boundary_centrality_score, -boundary_centrality_index, -boundary_centrality_ratio)
watershed_statistics <- watershed_statistics %>% select(-normalized_combined_boundary_length, -boundary_centrality_score, -boundary_centrality_index, -boundary_centrality_ratio)

# Export Spatial Watershed Dividedness variables
watershed_dividedness <- watersheds %>% left_join(results_normalized, by = "tnmid")
st_write(watershed_dividedness, "Data/watershed_dividedness.shp")

# Export spatial watershed statistics
watershed_stats_shp <- st_sf(watershed_statistics, geometry = st_geometry(watersheds))
st_write(watershed_stats_shp, "Data/watershed_statistics_full.shp")

# Remove coastal and wetland watersheds before export
watershed_statistics_shp <- watersheds %>% left_join(watershed_statistics, by = "tnmid")
watershed_statistics_shp <- watershed_statistics_shp %>% 
  filter(!watershed %in% watersheds_to_remove) %>%
  filter(!watershed == "Okefenokee Swamp")
st_write(watershed_statistics_shp, "Data/watershed_statistics.shp")

watershed_model_data <- watershed_statistics %>% 
  select(FLOOD_mean, intersection_count, boundary_centrality, dominance_ratio, 
         area_km2, watershed, pop_mean, AWATER_mean, severity_per_km2, 
         DROUGHT_mean, mean_elevation, DWATER_mean, wastewater_discharge_mean, watershed_dividedness) %>%
  filter(!watershed %in% watersheds_to_remove) %>%
  filter(!watershed == "Okefenokee Swamp")

# Begin regression
FloodModel <- lm(FLOOD_mean ~ intersection_count + 
                   boundary_centrality + 
                   pop_mean +
                   AWATER_mean +
                   # mean_elevation +
                   # severity_per_km2 +
                   # DROUGHT_mean +
                   log1p(wastewater_discharge_mean) +
                   # dominance_ratio + 
                   area_km2, data = watershed_model_data)

summary(FloodModel)

vars_used <- all.vars(formula(FloodModel))
na_rows <- watershed_statistics[!complete.cases(watershed_statistics[, vars_used]), ]

WaterModel <- lm(wastewater_discharge_mean ~ intersection_count + 
                   boundary_centrality + 
                   pop_mean +
                   AWATER_mean +
                   mean_elevation +
                   severity_per_km2 +
                   DROUGHT_mean +
                   dominance_ratio + 
                   area_km2, data = watershed_model_data)

summary(WaterModel)

# Begin regression
FireModel <- lm(severity_per_km2 ~ intersection_count + 
                  boundary_centrality + 
                  pop_mean +
                  AWATER_mean +
                  mean_elevation +
                  DROUGHT_mean +
                  dominance_ratio + 
                  area_km2, data = watershed_model_data)

summary(FireModel)

# Begin regression
DroughtModel <- lm(DROUGHT_mean ~ intersection_count + 
                     boundary_centrality + 
                     pop_mean +
                     AWATER_mean +
                     mean_elevation +
                     severity_per_km2 +
                     dominance_ratio + 
                     area_km2, data = watershed_model_data)

summary(DroughtModel)

# Phase 4: run tests on regression
# Multicolinearity
vif(FloodModel)

# Heteroskedasticity
bptest(FloodModel)

plot(FloodModel$fitted.values, resid(FloodModel),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
# Add linear trendline through residuals
trend <- lm(resid(FloodModel) ~ FloodModel$fitted.values)
abline(trend, col = "blue", lwd = 2)

nrow(watershed_statistics)
nrow(model.frame(FloodModel))

coeftest(FloodModel, vcov = vcovHC(FloodModel, type = "HC1"))

# Create a dataframe of fitted values and residuals
diagnostics <- data.frame(
  name = watershed_model_data$watershed,  # or whatever your watershed ID field is
  fitted = FloodModel$fitted.values,
  residual = resid(FloodModel)
)

extreme_point <- diagnostics %>%
  filter(abs(residuals) == max(abs(residuals)))  # or set a threshold like > 100

# View extreme values
head(diagnostics[order(diagnostics$residual), ], 5)  # most negative residuals
head(diagnostics[order(-diagnostics$residual), ], 5) # most positive residuals

ggplot(diagnostics, aes(x = fitted, y = residual)) +
  geom_point(color = "darkgray") +
  geom_smooth(method = "lm", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, color = "red") +
  theme_minimal() +
  labs(title = "Residuals vs Fitted with Trendline",
       x = "Fitted Values",
       y = "Residuals")

# Outliers
plot(cooks.distance(FloodModel))

ncvTest(FloodModel)  # Non-constant variance score test (alternative to bptest)

# Residuals vs. each predictor
par(mfrow = c(2, 3))
for (v in c("intersection_count", "boundary_centrality", "pop_mean", "AWATER_mean", "wastewater_discharge_mean", "area_km2", "mean_elevation")) {
  plot(watershed_model_data[[v]], resid(FloodModel), main = v,
       xlab = v, ylab = "Residuals")
  abline(h = 0, col = "red")
  print(v)
  print(skewness(watershed_model_data[[v]], na.rm = TRUE))
}
par(mfrow = c(1, 1))

# GWR
Flood_bw <- bw.gwr(FLOOD_mean ~ intersection_count + 
                     boundary_centrality + 
                     pop_mean +
                     AWATER_mean # + 
                     # mean_elevation +
                     # severity_per_km2 +
                     # DROUGHT_mean +
                      #area_km2 
                     , data = watershed_statistics_shp,
             adaptive = T) 
Flood_m <- gwr.basic(FLOOD_mean ~ intersection_count + 
                       boundary_centrality + 
                       pop_mean +
                       AWATER_mean, #+
                       # mean_elevation +
                       # severity_per_km2 +
                       # DROUGHT_mean +
                     #area_km2,
                        data = watershed_statistics_shp,
               adaptive = T,
               bw = Flood_bw)  

Flood_m
summary(Flood_m)
Flood_m$SDF
summary(Flood_m$SDF)
summary(Flood_m$lm)

AIC(Flood_m$lm)
vif(Flood_m$lm)
#WIF_res <- gwr.wif(Flood_m)
#summary(WIF_res)

# Extract fitted and residuals
residuals <- Flood_m$SDF$gwr.e
summary(residuals)
fitted_vals <- Flood_m$SDF$pred

plot(Flood_m$lm$fitted.values, resid(Flood_m$lm),
     xlab = "Fitted values",
     ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, col = "red")
# Add linear trendline through residuals
trend <- lm(resid(Flood_m$lm) ~ Flood_m$lm$fitted.values)
abline(trend, col = "blue", lwd = 2)

bptest(Flood_m$lm)

# plot
plot(fitted_vals, residuals,
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs Fitted (GWR)")
abline(h = 0, col = "red")

plot_vars <- c("intersection_count", "boundary_centrality", "pop_mean", "AWATER_mean")

Flood_long <- Flood_m$SDF %>%
  st_as_sf() %>%
  tidyr::pivot_longer(cols = all_of(plot_vars), 
                      names_to = "Variable", values_to = "Value")

st_write(Flood_long, "Data/Flood_long.shp")
st_write(Flood_m$SDF, "Data/Flood_SDF.shp")


tmap_mode('plot')

# Plot coefficients
fig = tm_shape(Flood_long) +
  tm_polygons(fill = "Value", 
              #size = 0.1, midpoint=NA,  
              fill.scale = tm_scale_intervals(
                style = "quantile",
                n = 4,
                values = c('#7D3C09', '#BF812E', '#35978F', '#003C30')
              ),
              fill.legend = tm_legend(title=c('Estimates on coeff for flood risk'),
              col.alpha = 0
  )) +
  tm_facets(by = "Variable", ncol = 2, nrow = 2)# +
  #tm_shape(watershed_statistics_shp) +  
  #tm_polygons(fill_alpha=0) + 
  tm_layout(
    legend.title.size = 1,
    legend.text.size = 0.7, 
    legend.position = c(0.02, 0.02),
    ) +
  tm_title("GWR local coeffs", size = 0.8, position = c(0.02, 0.27))

  summary(Flood_m$SDF$Local_R2)
  
fig

for (var in plot_vars) {
  map_data <- Flood_long %>% filter(Variable == var)
  
  fig <- tm_shape(map_data) +
    tm_polygons(
      fill = "Value",
      fill.scale = tm_scale_intervals(
        style = "quantile",
        n = 4,
        values = c('#7D3C09', '#BF812E', '#35978F', '#003C30')
      ),
      fill.legend = tm_legend(title = "Legend")
    ) +
    tm_layout(
      legend.title.size = 1,
      legend.text.size = 0.7,
      main.title = paste("Coefficient estimates for ", var),
      main.title.size = 1.0,
      legend.position = c("right", "top")
    )
  
  # Save each map to file
  tmap_save(fig, filename = paste0("Project/gwr_coeff_", var, ".png"), width = 6, height = 6, dpi = 300)
}

fig2 = tm_shape(Flood_m$SDF) + 
  tm_polygons(fill = "Local_R2") #+ #Local_R2
  #tm_shape(watershed_statistics_shp) + 
  #tm_polygons(col.alpha = 0)

fig2

fig3 = tm_shape(Flood_m$SDF) + 
  tm_polygons(fill = "residual") #+ #Local_R2
  #tm_shape(watershed_statistics_shp) + 
  #tm_polygons(col.alpha = 0)

fig3

#-Literature review - stakeholder workgroups used in divided watersheds and Army Corps water manuals for divided watersheds. 