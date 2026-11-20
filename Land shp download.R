install.packages("rnaturalearth")

library(sf)
library(rnaturalearth)
library(ggplot2)
library(dplyr)

land_world <- ne_download(
  scale = 10,
  type = "land",
  category = "physical",
  returnclass = "sf"
)

# Define study area
study_bbox <- st_bbox(
  c(
    xmin = 117,
    xmax = 121,
    ymin = 22,
    ymax = 26
  ),
  crs = st_crs(4326)
)

# Crop land polygons
land_TS <- st_crop(
  land_world,
  study_bbox
)

#check downloaded data
ggplot() +
  geom_sf(
    data = land_TS,
    fill = "grey70",
    color = "black",
    linewidth = 0.3
  ) +
  coord_sf(
    xlim = c(117, 121),
    ylim = c(22, 26),
    expand = FALSE
  ) +
  theme_bw()

#save as shp file
dir.create("land_Shapefiles", showWarnings = FALSE)

st_write(
  land_TS,
  "land_Shapefiles/Taiwan_Strait_Land.shp",
  delete_layer = TRUE
)
