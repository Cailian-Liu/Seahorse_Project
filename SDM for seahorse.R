#most codes are AI generated. I gave prompts and checked the codes.

library(terra)
library(ggplot2)
library(tidyverse)
library(tidyr)
library(readxl)
library(dplyr)
library(sf)
library(lubridate)


seahorse_presence_only <- read.csv("seahorse_presence_only.csv") 
#presence only data, for SDM
seahorse_pa_env_3month<- read.csv("seahorse_pa_env_3month.csv")
#data including all hauls (with/without seahorses), not all env cells
#for generating biased background points
seahorse_pa_env_3month %>%
  group_by(YearMonth) %>%
  summarise(n_cells=n()) #check How many sampled cells each month

env_SDM_selected <- read.csv("env_SDM_selected.csv")
#all env data
land_TS <- st_read("land_Shapefiles/Taiwan_Strait_Land.shp")
#land shp file


#generate biased background points
#use n_hauls in each cell to generate a sampling bias surface (KDE) with weighted cells, 
#and then get random points from the surface
#i.e. regions with higher n_haul should have more background points generated

#However, in SDM model, we have to use change abundance into presence-only data,
#Then using n_haul might overcorrect sampling bias. 
# using lg(n_haul + 1) to make the weight for each cell less extreme


#Plan A: no sampling bias correction.

#Plan B: using n_hauls as weight directly. 0 background points at no fishing region

#Plan C: using lg(n_haul + 1) as weight.

#Choose the best one as the weight in KDE and use weighted KDE to smooth sampling effort. 


# 1. generate background points----
fishing_effort_3month <- seahorse_pa_env_3month %>%
  select(
    YearMonth,
    cell_id,
    n_hauls)
background_pool <- env_SDM_selected %>%
  filter(
    !is.na(bottom_temp),
    !is.na(bottom_sal),
    !is.na(Depth)
  ) %>%
  left_join(
    fishing_effort_3month,
    by = c("YearMonth", "cell_id"),
    relationship = "one-to-one"
  ) %>%
  mutate(
    n_hauls = tidyr::replace_na(n_hauls, 0L)
  )
#check
background_pool %>%
  group_by(YearMonth) %>%
  summarise(
    total_cells = n(),
    surveyed_cells = sum(n_hauls > 0),
    unsurveyed_cells = sum(n_hauls == 0),
    total_hauls = sum(n_hauls),
    survey_coverage = surveyed_cells / total_cells * 100,
    .groups = "drop"
  )

# 4 different weights for each cell
background_pool <- background_pool %>%
  mutate(
    # M0: No sampling bias correction
    weight_uniform = 1,
    
    # M1: Raw haul count
    weight_raw = n_hauls,
    
    # M2: Square root transformation
    weight_sqrt = sqrt(n_hauls),
    
    # M3: Log transformation
    weight_log = log10(n_hauls + 1)
  ) %>%
  group_by(YearMonth) %>%
  mutate(
    prob_uniform = weight_uniform / sum(weight_uniform),
    prob_raw = weight_raw / sum(weight_raw),
    prob_sqrt = weight_sqrt / sum(weight_sqrt),
    prob_log = weight_log / sum(weight_log)
  ) %>%
  ungroup()


#uniform background points
n_background <- 500
set.seed(123)
background_uniform <- background_pool %>%
  group_by(YearMonth) %>%
  slice_sample(
    n = n_background,
    weight_by = prob_uniform,
    replace = TRUE #有放回采样
  ) %>%
  ungroup() %>%
  mutate(
    background_id = row_number(),
    Background_method = "Uniform",
    Response = 0L
  )

#raw n_haul as weight
set.seed(123)
background_raw <- background_pool %>%
  group_by(YearMonth) %>%
  slice_sample(
    n = n_background,
    weight_by = prob_raw,
    replace = TRUE
  ) %>%
  ungroup() %>%
  mutate(
    background_id = row_number(),
    Background_method = "Raw",
    Response = 0L
  )

#sqrt(n_haul) as weight
set.seed(123)
background_sqrt <- background_pool %>%
  group_by(YearMonth) %>%
  slice_sample(
    n = n_background,
    weight_by = prob_sqrt,
    replace = TRUE
  ) %>%
  ungroup() %>%
  mutate(
    background_id = row_number(),
    Background_method = "Sqrt",
    Response = 0L
  )

#lg(n_haul+1) as weight
set.seed(123)
background_log <- background_pool %>%
  group_by(YearMonth) %>%
  slice_sample(
    n = n_background,
    weight_by = prob_log,
    replace = TRUE
  ) %>%
  ungroup() %>%
  mutate(
    background_id = row_number(),
    Background_method = "Log",
    Response = 0L
  )

# combine all background points
background_all <- bind_rows(
  background_uniform,
  background_raw,
  background_sqrt,
  background_log
) %>%
  mutate(
    background_id = row_number(),
    Background_method = factor(
      Background_method,
      levels = c("Uniform", "Raw", "Sqrt", "Log")
    )
  )
#check how many cells were selected for background
background_all %>%
  group_by(YearMonth, Background_method) %>%
  summarise(
    total_background = n(),
    unique_cells = n_distinct(cell_id),
    repeated_background = n() - n_distinct(cell_id),
    background_in_surveyed = sum(n_hauls > 0),
    background_in_unsurveyed = sum(n_hauls == 0),
    .groups = "drop"
  ) %>%
  print(n = Inf)

#1.2 mapping 4 different background points----
background_count <- background_all %>%
  count(
    YearMonth,
    Background_method,
    cell_id,
    name = "n_background"  ) #count background points in each cell
background_map_all <- background_pool %>%
  crossing(
    Background_method = factor(
      c("Uniform", "Raw", "Sqrt", "Log"),
      levels = c("Uniform", "Raw", "Sqrt", "Log")
    )
  ) %>%
  left_join(
    background_count,
    by = c(
      "YearMonth",
      "Background_method",
      "cell_id"
    )
  ) %>%
  mutate(
    n_background = replace_na(n_background, 0L)
  ) #add unsampled cells into background cells

background_map_all %>%
  ggplot(
    aes(
      x = Longitude,
      y = Latitude,
      fill = n_background
    )
  ) +
  geom_tile(
    width = res_deg,
    height = res_deg
  ) +
  
  geom_sf(
    data = land_TS,
    fill = "grey65",
    color = "black",
    linewidth = 0.2,
    inherit.aes = FALSE
  ) +
  scale_fill_gradientn(
    colours = c(
      "white", "yellow", "gold",
      "orange", "red", "darkred"
    ),
    trans = "sqrt",
    breaks = c(0, 1, 5, 10, 20, 30, 50)
  ) +
  facet_grid(
    Background_method ~ YearMonth
  ) +
  
  scale_x_continuous(
    breaks = seq(117, 121, by = 2),
    labels = as.character
  ) +
  
  scale_y_continuous(
    breaks = seq(22, 26, by = 2),
    labels = as.character
  ) +
  
  coord_sf(
    crs = st_crs(4326),
    xlim = c(117, 121),
    ylim = c(22, 26),
    expand = FALSE
  ) +
  
  labs(
    x = "Longitude (°E)",
    y = "Latitude (°N)",
    fill = "Number of\nbackground points"
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold")
  )





#2. SDM build

#需要构建四套SDM数据，用于MaxEnt，GLM和randomforest
