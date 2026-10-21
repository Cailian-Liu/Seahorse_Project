library(terra)
library(ggplot2)
library(readxl)
library(dplyr)
library(sf)


bottom_env_all <- read.csv("./Copernicus_TaiwanStrait/near_bottom/bottom_env_all.csv")
logbook_data <- read_excel("./Seahorse data/Logbook data.xlsx")
str(logbook_data)
logbook_data <- logbook_data %>%
  mutate(
    Date = as.Date(paste(Month, Day_corrected, sep = "-"))
  )
sum(is.na(logbook_data$Date))

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






