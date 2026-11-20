#most codes are AI generated. I gave prompts and checked the codes.

library(terra)
library(ggplot2)
library(tidyverse)
library(readxl)
library(dplyr)
library(sf)
library(lubridate)


bottom_env_monthly <- read.csv("./Copernicus_TaiwanStrait/near_bottom/bottom_env_monthly.csv")
str(bottom_env_monthly)
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

land_TS <- st_read(
  "land_Shapefiles/Taiwan_Strait_Land.shp")


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
depth <- depth %>%
  filter(
    Longitude >= 117,
    Longitude <= 121,
    Latitude >= 22,
    Latitude <= 26  ) 

ggplot(depth, aes(Longitude, Latitude, fill = Depth2)) +
  geom_tile(
    width = 5 / 60,
    height = 5 / 60
  ) +
  scale_fill_viridis_c(name = "Depth") +
  coord_equal() +
  theme_bw() #是规则的网格数据


#combine depth + bottom env data
# 设置显示精度
options(digits = 12)
# Copernicus数据的经纬度
head(sort(unique(bottom_env_monthly$Longitude)), 10)
head(sort(unique(bottom_env_monthly$Latitude)), 10)
# Depth数据的经纬度
head(sort(unique(depth$Longitude)), 10)
head(sort(unique(depth$Latitude)), 10)

# 检查两个数据集的空间范围
range(bottom_env_monthly$Longitude)
range(bottom_env_monthly$Latitude)
range(depth$Longitude)
range(depth$Latitude)

#check if the data are in the same grid
res_deg <- 5 / 60
# 检查坐标与标准5′网格中心的最大偏差
grid_error <- function(x) {
  max(abs(x - round(x / res_deg) * res_deg),
      na.rm = TRUE)
}
grid_error(bottom_env_monthly$Longitude)
grid_error(bottom_env_monthly$Latitude)
grid_error(depth$Longitude)
grid_error(depth$Latitude)
#all error are very small, good to keep going

#创建统一的空间索引，用于后续数据合并
# Copernicus环境数据
env_monthly_grid <- bottom_env_monthly %>%
  mutate(
    lon_id = round(Longitude / res_deg),
    lat_id = round(Latitude / res_deg)
  )

# 水深数据
depth_grid <- depth %>%
  mutate(
    lon_id = round(Longitude / res_deg),
    lat_id = round(Latitude / res_deg)
  ) %>%
  select(lon_id, lat_id, Depth = Depth1)
#combine 2 data
env_monthly_all <- env_monthly_grid %>%
  left_join(
    depth_grid,
    by = c("lon_id", "lat_id"),
    relationship = "many-to-one"
  ) %>%
  select(-lon_id, -lat_id)
str(env_monthly_all)

# 检查合并前后行数是否一致
nrow(bottom_env_monthly)
nrow(env_monthly_all)
# 检查有多少记录匹配到了水深
env_monthly_all %>%
  summarise(
    total_records = n(),
    matched = sum(!is.na(Depth)),
    unmatched = sum(is.na(Depth)),
    matched_percentage = mean(!is.na(Depth)) * 100
  )

#plot missing cells:all lands
env_monthly_all %>%
  filter(YearMonth == "2022-08") %>%
  ggplot(aes(x = Longitude,
             y = Latitude,
             fill = is.na(Depth))) +
  geom_tile(width = res_deg, height = res_deg) +
  scale_fill_manual(
    values = c("FALSE" = "steelblue", "TRUE" = "grey90"),
    labels = c("Matched", "Missing"),
    name = "Depth data"
  ) +
  coord_quickmap() +
  theme_bw()

#final env data for SDM: keep only data from ocean
env_SDM <- env_monthly_all %>%
  filter(
    !is.na(Depth),
    !is.na(bottom_temp),
    !is.na(bottom_sal)
  ) %>%
  select(
    YearMonth,
    Longitude,
    Latitude,
    bottom_temp,
    bottom_sal,
    Depth
  )
str(env_SDM)

# 检查每个变量是否还有NA
colSums(is.na(env_SDM))
write.csv(env_SDM,"env_SDM.csv", row.names = FALSE)



#change seahorse data to presence/absence data, according to grids from env_SDM

#create grids based on Copernicus Marine data
# Get the coordinates of grid centers
lon_min <- min(bottom_env_monthly$Longitude)
lon_max <- max(bottom_env_monthly$Longitude)
lat_min <- min(bottom_env_monthly$Latitude)
lat_max <- max(bottom_env_monthly$Latitude)
# Create reference raster
ref_r <- rast(
  xmin = lon_min - res_deg / 2,
  xmax = lon_max + res_deg / 2,
  ymin = lat_min - res_deg / 2,
  ymax = lat_max + res_deg / 2,
  resolution = res_deg,
  crs = "EPSG:4326"
)
ref_r
#create cell_id
env_SDM <- env_SDM %>%
  mutate(
    cell_id = cellFromXY(
      ref_r,
      cbind(Longitude, Latitude)
    )
  )
logbook_grid <- logbook_data %>%
  mutate(
    haul_id = row_number(),
    YearMonth = format(Date, "%Y-%m"),
    cell_id = cellFromXY(
      ref_r,
      cbind(Longitude_haulback, Latitude_haulback) 
      #use location when hauling back as seahorse location
    )
  )

#generate monthly presence/absence data
seahorse_pa_monthly <- logbook_grid %>%
  filter(!is.na(cell_id)) %>%
  group_by(YearMonth, cell_id) %>%
  summarise(
    n_hauls = n(), # #of hauls in each cell, useful for generating biased background points
    
    n_positive_hauls = sum(
      Seahorse > 0, na.rm = TRUE # #of hauls with seahorse, useless
    ),
    
    total_duration = sum(
      Duration, na.rm = TRUE # fishing duration in this cell, useless
    ),
    
    total_seahorse = if (all(is.na(Seahorse))) {
      NA_real_
    } else {
      sum(Seahorse, na.rm = TRUE) # # of seahorses in each cell, useless
    },
    
    Presence = case_when( #presence of seahorse in each cell, useful
      any(Seahorse > 0, na.rm = TRUE) ~ 1L,
      all(!is.na(Seahorse) & Seahorse == 0) ~ 0L,
      TRUE ~ NA_integer_
    ),
    
    .groups = "drop"
  )
str(seahorse_pa_monthly)


#combine seahorse data and env data
seahorse_pa_monthly_env <- seahorse_pa_monthly %>%
  left_join(
    env_SDM,
    by = c("YearMonth", "cell_id"),
    relationship = "many-to-one" #
  ) %>%
  mutate(
    valid_env = complete.cases(
      bottom_temp,
      bottom_sal,
      Depth
    )
  )

#sample size calculation
monthly_sample_size <- seahorse_pa_monthly_env %>%
  group_by(YearMonth) %>%
  summarise(
    # Total number of hauls
    total_hauls = sum(n_hauls),
    
    # Presence and absence before environmental filtering
    presence_cells = sum(Presence == 1, na.rm = TRUE),
    absence_cells = sum(Presence == 0, na.rm = TRUE),
    
    # Cells with complete environmental data
    valid_surveyed_cells = sum(
      valid_env & !is.na(Presence)
    ),
    
    valid_presence_cells = sum(
      Presence == 1 & valid_env, na.rm = TRUE
    ),
    
    valid_absence_cells = sum(
      Presence == 0 & valid_env, na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    presence_percentage = if_else(
      valid_surveyed_cells > 0,
      valid_presence_cells / valid_surveyed_cells * 100,
      NA_real_
    ),
    
    enough_samples = valid_presence_cells > 30 
    #only August, September, and October have >30 sampled grids
  )

print(monthly_sample_size, n = Inf, width = Inf)
#only August, September, and October have >30 sampled grids

#bar plot of sample size 
monthly_sample_size %>%
  select(
    YearMonth,
    valid_presence_cells,
    valid_absence_cells
  ) %>%
  pivot_longer(
    cols = -YearMonth,
    names_to = "Type",
    values_to = "n_cells"
  ) %>%
  mutate(
    Type = factor(
      Type,
      levels = c(
        "valid_absence_cells",
        "valid_presence_cells"
      ),
      labels = c("Absence", "Presence")
    )
  ) %>%
  ggplot(
    aes(x = YearMonth, y = n_cells, fill = Type)
  ) +
  geom_col(width = 0.75) +
  geom_hline(
    yintercept = 30,
    linetype = "dashed",
    color = "red"
  ) +
  scale_fill_manual(
    values = c(
      "Presence" = "steelblue",
      "Absence" = "grey75"
    )
  ) +
  labs(
    x = "Month",
    y = "Number of surveyed cells",
    fill = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(
      angle = 45, hjust = 1
    )
  )

#mapping presence/absence/unsampled cells
seahorse_SDM_full <- env_SDM %>%
  left_join(
    seahorse_pa_monthly,
    by = c("YearMonth", "cell_id")
  ) %>%
  mutate(
    Status = case_when(
      Presence == 1 ~ "Presence",
      Presence == 0 ~ "Absence",
      is.na(n_hauls) ~ "Unsampled",
      TRUE ~ "Unknown"
    ),
    n_hauls = replace_na(n_hauls, 0L)
  )
seahorse_SDM_full %>%
  filter(
    YearMonth %in% seahorse_pa_monthly$YearMonth
  ) %>%
  ggplot(
    aes(x = Longitude, y = Latitude, fill = Status)
  ) +
  geom_tile(
    width = res_deg,
    height = res_deg
  ) +
  # Land and coastline
  geom_sf(
    data = land_TS,
    fill = "grey65",
    color = "black",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(
    values = c(
      "Presence" = "steelblue",
      "Absence" = "orange",
      "Unsampled" = "grey90",
      "Unknown" = "purple"
    )
  ) +
  facet_wrap(~ YearMonth, ncol = 4) +
  labs(
    x = "Longitude",
    y = "Latitude",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(
      fill = "white"
    ),
    legend.position = "bottom",
    strip.background = element_rect(
      fill = "grey90"
    ),
    strip.text = element_text(
      face = "bold"
    )
  )

#keep data from only three month with sample>30
selected_months <- c("2022-08", "2022-09", "2022-10")
seahorse_pa_env_3month <- seahorse_pa_monthly_env %>%
  filter(
    YearMonth %in% selected_months,
    !is.na(Presence),
    valid_env  )
write.csv(seahorse_pa_env_3month,"seahorse_pa_env_3month.csv", row.names = FALSE)

#check # of hauls (with/without seahorse) in each cell for this whole region
seahorse_pa_env_3month %>%
  filter(n_hauls > 0) %>%
  group_by(YearMonth) %>%
  summarise(
    n_cells = n(),
    total_hauls = sum(n_hauls),
    min_hauls = min(n_hauls),
    max_hauls = max(n_hauls),
    mean_hauls = mean(n_hauls),
    median_hauls = median(n_hauls),
    sd_hauls = sd(n_hauls),
    .groups = "drop"
  )
#histgram of # of hauls in each cell
seahorse_pa_env_3month %>%
  ggplot(aes(x = n_hauls)) +
  geom_histogram(
    binwidth = 1,
    center = 1,
    fill = "steelblue",
    color = "black"
  ) +
  facet_wrap(~ YearMonth, ncol = 3) +
  labs(
    x = "Number of hauls per cell",
    y = "Number of cells"
  ) +
  theme_bw(base_size = 14)

seahorse_pa_env_3month %>%
  filter(n_hauls > 0) %>%
  ggplot(aes(x = YearMonth, y = n_hauls)) +
  geom_boxplot(fill = "steelblue") +
  geom_jitter(
    width = 0.15,
    alpha = 0.3,
    size = 1.5
  ) +
  labs(
    x = "Month",
    y = "Number of hauls per cell"
  ) +
  theme_bw(base_size = 14)

#only hauls with seahorse
seahorse_pa_env_3month %>%
  filter(n_positive_hauls > 0) %>%
  group_by(YearMonth) %>%
  summarise(
    n_cells = n(),
    total_positive_hauls = sum(n_positive_hauls),
    min_positive_hauls = min(n_positive_hauls),
    max_positive_hauls = max(n_positive_hauls),
    mean_positive_hauls = mean(n_positive_hauls),
    median_positive_hauls = median(n_positive_hauls),
    sd_positive_hauls = sd(n_positive_hauls),
    .groups = "drop"
  ) %>%
  print(width = Inf)

seahorse_pa_env_3month %>%
  ggplot(aes(x = n_positive_hauls)) +
  geom_histogram(
    binwidth = 1,
    center = 1,
    fill = "steelblue",
    color = "black"
  ) +
  facet_wrap(~ YearMonth, ncol = 3) +
  labs(
    x = "Number of positive hauls per cell",
    y = "Number of cells"
  ) +
  theme_bw(base_size = 14)


#other dataset needed for following SDM
seahorse_presence_only <- seahorse_pa_env_3month %>%
  filter(Presence == 1)
write.csv(seahorse_presence_only,"seahorse_presence_only.csv", row.names = FALSE)

env_SDM_selected <- env_SDM %>%
  filter(YearMonth %in% selected_months)
write.csv(env_SDM_selected, "env_SDM_selected.csv", row.names = FALSE)
