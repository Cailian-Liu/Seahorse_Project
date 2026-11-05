library(terra)
library(ggplot2)
library(readxl)
library(dplyr)
library(sf)
library(lubridate)


bottom_env_all <- read.csv("./Copernicus_TaiwanStrait/near_bottom/bottom_env_all.csv")
logbook_data <- read_excel("./Seahorse data/Logbook data.xlsx")
str(logbook_data)
logbook_data <- logbook_data %>%
  mutate(
    Date = as.Date(paste(Month, Day_corrected, sep = "-"))
  )
range(logbook_data$Date)
sum(is.na(logbook_data$Date))
sum(is.na(logbook_data$Duration))
sum(logbook_data$Seahorse != 0, na.rm = TRUE)


#### Depth data
depth <- read.table(
  "cs_depth.txt",
  header = FALSE,
  sep = "\t")
colnames(depth) <- c(
  "Longitude",
  "Latitude",
  "Depth1",
  "Depth2")
head(depth)
ggplot(depth, aes(Longitude, Latitude, fill = Depth2)) +
  geom_tile(
    width = 5 / 60,
    height = 5 / 60
  ) +
  scale_fill_viridis_c(name = "Depth") +
  coord_equal() +
  theme_bw() #是规则的网格数据

#convert to raster
#恢复规则网格数据
res_deg <- 5 / 60
depth_fixed <- depth %>%
  mutate(
    Longitude = round(Longitude / res_deg) * res_deg,
    Latitude  = round(Latitude / res_deg) * res_deg  )
max(abs(depth$Longitude - depth_fixed$Longitude)) #坐标修正幅度小
max(abs(depth$Latitude - depth_fixed$Latitude))
depth_fixed %>%
  count(Longitude, Latitude) %>%
  filter(n > 1) #无重复坐标

depth_r <- rast(
  depth_fixed[, c("Longitude", "Latitude", "Depth2")],
  type = "xyz",
  crs = "EPSG:4326") #create raster
print(depth_r)
res(depth_r) # 5 arc-minutes


####Near-bottom environment data###
#also 5 arc-minutes resolution

#logbook data resolution: irregular spatial points


#将depth data跟env data网格对齐--这只是针对一天的
env_day <- bottom_env_all %>%
  filter(Date == "2022-08-01")
#手动创建5‘raster
res_deg <- 5 / 60
lon0 <- min(env_day$Longitude, na.rm = TRUE)
lat0 <- min(env_day$Latitude, na.rm = TRUE)
# Grid indices, starting at zero
col_id <- round((env_day$Longitude - lon0) / res_deg)
row_id <- round((env_day$Latitude - lat0) / res_deg)
# Check whether the coordinates fit a 5' grid --误差都很小，可以
max(abs(env_day$Longitude - (lon0 + col_id * res_deg)))
max(abs(env_day$Latitude - (lat0 + row_id * res_deg)))
env_template <- rast(
  ncols = max(col_id) + 1,
  nrows = max(row_id) + 1,
  xmin = lon0 - res_deg / 2,
  xmax = lon0 + (max(col_id) + 0.5) * res_deg,
  ymin = lat0 - res_deg / 2,
  ymax = lat0 + (max(row_id) + 0.5) * res_deg,
  crs = "EPSG:4326")
res(env_template)
ext(env_template)
# Find raster cell for each original observation
cell_id <- cellFromXY(
  env_template,
  cbind(
    lon0 + col_id * res_deg,
    lat0 + row_id * res_deg
  )
)
# Check for duplicate cells
sum(duplicated(cell_id))

#将所有near-bottom环境变量转换成raster
# Near-bottom environmental variables
env_vars <- c(
  "bottom_depth",
  "bottom_temp",
  "bottom_sal",
  "bottom_u",
  "bottom_v",
  "bottom_speed")
# Create one raster layer for each variable
env_layers <- lapply(env_vars, function(v) {
  r <- env_template
  values(r) <- NA_real_
  r[cell_id] <- env_day[[v]]
  names(r) <- v
  return(r)
})
# Combine layers
env_r <- rast(env_layers)
names(env_r)
nlyr(env_r)

#将depth_r对齐到Copernicus网格
depth_aligned <- resample(
  depth_r,
  env_template,
  method = "bilinear")
names(depth_aligned) <- "depth_external"
compareGeom(
  depth_aligned,
  env_r,
  stopOnError = FALSE) #check
env_all_r <- c(
  env_r,
  depth_aligned)
names(env_all_r)
#save as GEOTIFF, multiple rasters
#writeRaster(env_all_r,"Copernicus_TaiwanStrait/near_bottom/env_all_2022-08-01.tif",  overwrite = TRUE)

#read 
#env_all_r <- rast("Copernicus_TaiwanStrait/near_bottom/env_all_2022-08-01.tif")
names(env_all_r)
res(env_all_r)
ext(env_all_r)
crs(env_all_r)

#check how many hauls in one 5 arc-minutes cell
haul_vect <- vect(
  logbook,
  geom = c(
    "Longitude_haulback",
    "Latitude_haulback"
  ),
  crs = "EPSG:4326"
)
# Assign each haul to a raster cell
logbook$env_cell <- cellFromXY(
  env_all_r,
  crds(haul_vect))
# Number of hauls per cell
haul_density <- logbook %>%
  count(env_cell, name = "n_hauls")
hist(haul_density$n_hauls)
haul_count_r <- rasterize(
  haul_vect,
  env_all_r[[1]],
  field = 1,
  fun = "sum",
  background = 0)
plot(haul_count_r) #整年的数据，确实有少数集中在某区域

#check how many hauls in one cell using only August data
logbook_env$cell_id <- cellFromXY(
  env_all_r,
  as.matrix(
    logbook_env[, c(
      "Longitude_haulback",
      "Latitude_haulback"
    )]
  )
)
# 1. Select hauls on August 1, 2022
haul_day <- logbook_env %>%
  filter(Month == ("2022-08"))
# 2. Match hauls to reference grid cells
haul_day$cell_id <- cellFromXY(
  env_all_r[[1]],
  as.matrix(
    haul_day[, c("Longitude_haulback", "Latitude_haulback")]
  )
)

# 3. Count hauls per cell
haul_count <- haul_day %>%
  filter(!is.na(cell_id)) %>%
  count(cell_id, name = "n_hauls")

# 4. Create raster
haul_r <- env_all_r[[1]]
values(haul_r) <- 0

haul_r[haul_count$cell_id] <- haul_count$n_hauls
names(haul_r) <- "n_hauls"

# 5. Plot
plot(
  haul_r,
  main = "Number of hauls per 5' cell (2022-08)",
  col = hcl.colors(20, "YlOrRd"),
  axes = TRUE
)

#




















###### Env + Seahorse data merge ----

#计算地理距离，找到每网最近的有有效环境数据的网格，再根据日期读取当天的环境值
# 1. Prepare logbook data
logbook <- logbook_data %>%
  mutate(
    haul_id = row_number(),
    Date = as.Date(Date)
  )

# 2. Prepare environmental data
env <- bottom_env_all %>%
  mutate(Date = as.Date(Date)) %>%
  filter(
    !is.na(bottom_depth),
    !is.na(bottom_temp),
    !is.na(bottom_sal),
    !is.na(bottom_u),
    !is.na(bottom_v)
  )

# 3. Convert to sf points (WGS84)
haul_sf <- st_as_sf(
  logbook %>%
    filter(
      !is.na(Longitude_haulback),
      !is.na(Latitude_haulback)
    ),
  coords = c("Longitude_haulback", "Latitude_haulback"),
  crs = 4326,
  remove = FALSE
)

# 4. Create unique environmental grid points
env_grid <- env %>%
  distinct(Longitude, Latitude)

env_sf <- st_as_sf(
  env_grid,
  coords = c("Longitude", "Latitude"),
  crs = 4326,
  remove = FALSE
)

# 5. Find nearest environmental grid for each haul
nearest_id <- st_nearest_feature(haul_sf, env_sf)

# Calculate distance to matched grid (meters)
distance_m <- as.numeric(
  st_distance(
    haul_sf,
    env_sf[nearest_id, ],
    by_element = TRUE
  )
)

# 6. Attach grid coordinates to each haul
haul_matched <- haul_sf %>%
  st_drop_geometry() %>%
  mutate(
    env_lon = env_grid$Longitude[nearest_id],
    env_lat = env_grid$Latitude[nearest_id],
    match_distance_m = distance_m
  )

# 7. Join daily environmental values
env_for_join <- env %>%
  rename(
    env_lon = Longitude,
    env_lat = Latitude
  )

logbook_env <- haul_matched %>%
  left_join(
    env_for_join,
    by = c("Date", "env_lon", "env_lat"),
    relationship = "many-to-one"
  )

# # 8. Restore hauls with missing coordinates：上一步后已经有了所有网次
# logbook_env <- logbook %>%
#   select(haul_id) %>%
#   left_join(
#     logbook_env,
#     by = "haul_id",
#     relationship = "one-to-one"
#   ) %>%
#   arrange(haul_id)

summary(logbook_env$match_distance_m)
#Copernicus Marine resolution is 1/12(0.0833) degree ~ 8-9km

# 10. Export
write.csv(
  logbook_env,
  "Seahorse data/log_env.csv",
  row.names = FALSE
)


##### Check how many hauls share exactly same env from same location&date ----
env_duplicates <- logbook_env %>%
  filter(
    !is.na(Date),
    !is.na(env_lon),
    !is.na(env_lat),
    !is.na(bottom_temp)
  ) %>%
  group_by(Date, env_lon, env_lat) %>%
  summarise(
    n_hauls = n(),
    .groups = "drop"
  )

# Overall summary
env_duplicates %>%
  summarise(
    n_env_groups = n(),
    n_unique_groups = sum(n_hauls == 1),
    n_repeated_groups = sum(n_hauls > 1),
    n_hauls_shared = sum(n_hauls[n_hauls > 1]),
    n_extra_duplicates = sum(n_hauls - 1),
    max_hauls_per_group = max(n_hauls)
  )
# 只有1网使用的环境数据：758个；被至少2网共享的环境组合数：440个
# 同一天同一网格最多匹配了多少网：10
env_duplicates %>%
  count(n_hauls, name = "n_env_groups") #查看每个环境组合被多少网共享





###### haul location distance check ----

#第一种距离：每网与其他网次的最近邻距离
# Convert haulback coordinates to sf
haul_sf <- logbook_data %>%
  mutate(haul_id = row_number()) %>%
  filter(
    !is.na(Longitude_haulback),
    !is.na(Latitude_haulback)
  ) %>%
  st_as_sf(
    coords = c("Longitude_haulback", "Latitude_haulback"),
    crs = 4326,
    remove = FALSE
  )
# Calculate pairwise geographic distances (meters)
dist_matrix <- st_distance(haul_sf)
# Exclude distance from a haul to itself
diag(dist_matrix) <- Inf
# Find nearest neighbor for each haul
nearest_distance <- apply(
  as.matrix(dist_matrix),
  1,
  min
)
# Store results in km
haul_sf$nearest_distance_km <- as.numeric(
  nearest_distance
) / 1000
# Summary
summary(haul_sf$nearest_distance_km) #0-19.8km, 平均=0.88km， median=0.54km
ggplot(
  st_drop_geometry(haul_sf),
  aes(x = nearest_distance_km)
) +
  geom_histogram(
    binwidth = 0.5,
    fill = "steelblue",
    color = "white"
  ) +
  labs(
    x = "Nearest-neighbor distance (km)",
    y = "Number of hauls",
    title = "Spatial Distance Between Seahorse Trawl Hauls"
  ) +
  theme_bw() #非常右偏

#第二种距离：同一艘船、同一天相邻网次之间的距离
#肯定非常近，没有算的必要。






