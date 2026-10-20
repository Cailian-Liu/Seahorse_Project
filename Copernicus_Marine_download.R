
# install.packages(c("reticulate"))

library(reticulate)
library(terra)
library(ggplot2)
library(reticulate)

# Use the Python environment that worked previously
cm <- import("copernicusmarine")

# Output folder
out_dir <- "Copernicus_TaiwanStrait"
dir.create(out_dir, showWarnings = FALSE,
           recursive = TRUE)

# Generate 13 monthly start dates
months <- seq(
  as.Date("2021-10-01"),
  as.Date("2022-10-31"),
  by = "month"
)

# Store download results
download_log <- data.frame(
  Month = character(),
  Status = character(),
  Message = character()
)

#download data for 13 months
for (i in seq_along(months)) {
  
  month_start <- months[i]
  
  month_end <- seq(
    month_start,
    by = "month",
    length.out = 2
  )[2] - 1
  
  # First and last day of each month
  month_end <- seq(
    month_start, by = "month", length.out = 2
  )[2] - 1
  
  month_label <- format(month_start, "%Y_%m")
  
  filename <- paste0(
    "ocean_", month_label, ".nc"
  )
  
  filepath <- file.path(out_dir, filename)
  
  # Skip previously downloaded files
  if (file.exists(filepath) &&
      file.size(filepath) > 0) {
    
    cat("Skipping existing:", filename, "\n")
    
    download_log <- rbind(
      download_log,
      data.frame(
        Month = month_label,
        Status = "Skipped",
        Message = "File already exists"
      )
    )
    
    next
  }
  
  cat("\nDownloading:", month_label, "\n")
  
  # Download with error handling
  result <- tryCatch({
    
    cm$subset(
      dataset_id =
        "cmems_mod_glo_phy_my_0.083deg_P1D-m",
      
      variables = list(
        "thetao", "so", "uo", "vo"
      ),
      
      minimum_longitude = 117,
      maximum_longitude = 119.5,
      minimum_latitude = 22.5,
      maximum_latitude = 24.5,
      
      start_datetime = paste0(
        month_start, "T00:00:00"
      ),
      
      end_datetime = paste0(
        month_end, "T23:59:59"
      ),
      
      output_directory = out_dir,
      output_filename = filename,
      
      overwrite = TRUE
    )
    
    "Success"
    
  }, error = function(e) {
    
    cat("Download failed:", conditionMessage(e), "\n")
    
    paste0("Failed: ", conditionMessage(e))
  })
  
  download_log <- rbind(
    download_log,
    data.frame(
      Month = month_label,
      Status = ifelse(
        result == "Success", "Success", "Failed"
      ),
      Message = result
    )
  )
}

# Display results: 检查下载的数据
print(download_log)

files <- sort(list.files(
  "Copernicus_TaiwanStrait",
  pattern = "^ocean_[0-9]{4}_[0-9]{2}\\.nc$",
  full.names = TRUE
))

file_check <- lapply(files, function(f) {
  
  nc <- nc_open(f)
  
  dates <- as.Date(
    as.POSIXct(
      ncvar_get(nc, "time") * 3600,
      origin = "1950-01-01",
      tz = "UTC"
    )
  )
  
  result <- data.frame(
    File = basename(f),
    Start = min(dates),
    End = max(dates),
    Days = length(dates),
    Depth_levels = nc$dim$depth$len,
    Variables_OK = all(
      c("thetao", "so", "uo", "vo") %in%
        names(nc$var)
    )
  )
  
  nc_close(nc)
  result
})
file_check <- do.call(rbind, file_check)
print(file_check)

###### get near bottom data #################################
# ==================================================
# 1. Set input/output directories
# ==================================================

data_dir <- "Copernicus_TaiwanStrait"

output_dir <- file.path(data_dir, "near_bottom")
dir.create(output_dir, showWarnings = FALSE,
           recursive = TRUE)

# Find all monthly NetCDF files
nc_files <- sort(list.files(
  data_dir,
  pattern = "^ocean_[0-9]{4}_[0-9]{2}\\.nc$",
  full.names = TRUE
))

print(basename(nc_files))
cat("Number of files:", length(nc_files), "\n")

stopifnot(length(nc_files) == 13)
# ==================================================
# 2. Function to process one monthly NetCDF
# ==================================================
extract_near_bottom <- function(nc_file) {
  
  cat("\nProcessing:", basename(nc_file), "\n")
  
  nc <- nc_open(nc_file)
  
  # Ensure NetCDF closes even if an error occurs
  on.exit(nc_close(nc), add = TRUE)
  
  # Read coordinates
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  depth <- ncvar_get(nc, "depth")
  time_hours <- ncvar_get(nc, "time")
  
  dates <- as.Date(
    as.POSIXct(
      time_hours * 3600,
      origin = "1950-01-01",
      tz = "UTC"
    )
  )
  
  # Read four 4D variables
  # ncdf4 applies scale_factor and add_offset
  # and converts _FillValue to NA
  temp <- ncvar_get(nc, "thetao",
                    collapse_degen = FALSE)
  sal <- ncvar_get(nc, "so",
                   collapse_degen = FALSE)
  u <- ncvar_get(nc, "uo",
                 collapse_degen = FALSE)
  v <- ncvar_get(nc, "vo",
                 collapse_degen = FALSE)
  
  # Create longitude-latitude-time grid
  grid <- expand.grid(
    lon_id = seq_along(lon),
    lat_id = seq_along(lat),
    time_id = seq_along(dates)
  )
  
  # Extract deepest valid layer
  get_bottom <- function(i, j, t) {
    
    valid <- which(
      is.finite(temp[i, j, , t]) &
        is.finite(sal[i, j, , t]) &
        is.finite(u[i, j, , t]) &
        is.finite(v[i, j, , t])
    )
    
    if (length(valid) == 0) {
      return(c(
        bottom_depth = NA_real_,
        bottom_temp = NA_real_,
        bottom_sal = NA_real_,
        bottom_u = NA_real_,
        bottom_v = NA_real_
      ))
    }
    
    # Deepest valid model layer
    k <- valid[which.max(depth[valid])]
    
    c(
      bottom_depth = depth[k],
      bottom_temp = temp[i, j, k, t],
      bottom_sal = sal[i, j, k, t],
      bottom_u = u[i, j, k, t],
      bottom_v = v[i, j, k, t]
    )
  }
  
  # Apply extraction across all cells and dates
  bottom_values <- t(
    mapply(
      get_bottom,
      grid$lon_id,
      grid$lat_id,
      grid$time_id
    )
  )
  
  # Combine coordinates and values
  bottom_env <- bind_cols(
    data.frame(
      Date = dates[grid$time_id],
      Longitude = lon[grid$lon_id],
      Latitude = lat[grid$lat_id]
    ),
    as.data.frame(bottom_values)
  ) %>%
    mutate(
      bottom_speed = sqrt(
        bottom_u^2 + bottom_v^2
      )
    )
  
  return(bottom_env)
}


# ==================================================
# 3. Process all 13 months
# ==================================================
monthly_results <- list()

for (i in seq_along(nc_files)) {
  
  nc_file <- nc_files[i]
  
  month_label <- sub(
    "^ocean_([0-9]{4}_[0-9]{2})\\.nc$",
    "\\1",
    basename(nc_file)
  )
  
  # Extract near-bottom data
  result <- extract_near_bottom(nc_file)
  
  # Save monthly CSV
  output_file <- file.path(
    output_dir,
    paste0("bottom_env_", month_label, ".csv")
  )
  
  write.csv(
    result,
    output_file,
    row.names = FALSE
  )
  
  monthly_results[[i]] <- result
  
  cat(
    "Finished:", month_label,
    "| Rows:", nrow(result), "\n"
  )
}


# ==================================================
# 4. Combine all months
# ==================================================

bottom_env_all <- bind_rows(monthly_results) %>%
  arrange(Date, Latitude, Longitude)

# Save combined dataset
write.csv(
  bottom_env_all,
  file.path(
    output_dir,
    "bottom_env_all.csv"
  ),
  row.names = FALSE
)

cat("\nAll months combined successfully!\n")

# ==================================================
# 5. Check results
# ==================================================

print(head(bottom_env_all))
print(summary(bottom_env_all))
#bottom speed: uo+vo，表示日平均水平流速大小
#bottom depth: 模型中最深有效水层的中心深度，不是实际海底水深

cat("Date range:\n")
print(range(bottom_env_all$Date))

cat("Number of days:",
    n_distinct(bottom_env_all$Date), "\n")

cat("Number of grid cells:",
    n_distinct(
      bottom_env_all$Longitude,
      bottom_env_all$Latitude
    ), "\n")


bottom_env_all %>%
  summarise(
    temp_min = min(bottom_temp, na.rm = TRUE),
    temp_max = max(bottom_temp, na.rm = TRUE),
    sal_min = min(bottom_sal, na.rm = TRUE),
    sal_max = max(bottom_sal, na.rm = TRUE),
    speed_max = max(bottom_speed, na.rm = TRUE),
    depth_max = max(bottom_depth, na.rm = TRUE)
  )

depth_check <- bottom_env_all %>%
  group_by(Longitude, Latitude) %>%
  summarise(
    n_depth = n_distinct(
      bottom_depth
    ),
    .groups = "drop"
  )

depth_check %>%
  count(n_depth) 

#检查396天中，有多少天存在有效数据
grid_missing <- bottom_env_all %>%
  group_by(Longitude, Latitude) %>%
  summarise(
    n_days = n_distinct(Date),
    n_valid = sum(!is.na(bottom_depth)),
    n_missing = sum(is.na(bottom_depth)),
    missing_prop = n_missing / n(),
    .groups = "drop"
  )

# Summarize missing values
grid_missing %>%
  count(n_missing)

#plot哪些网格缺失数据：缺失数据的是陆地
grid_missing <- grid_missing %>%
  mutate(
    Status = case_when(
      missing_prop == 0 ~ "Complete",
      missing_prop == 1 ~ "All missing",
      TRUE ~ "Partially missing"
    )
  )

ggplot(
  grid_missing,
  aes(x = Longitude, y = Latitude)
) +
  geom_tile(
    aes(fill = Status),
    width = 1/12,
    height = 1/12
  ) +
  scale_fill_manual(
    values = c(
      "Complete" = "#2A9D8F",
      "All missing" = "#E76F51",
      "Partially missing" = "#F4A261"
    ),
    name = "Data availability"
  ) +
  coord_fixed(
    ratio = 1 / cos(23.5 * pi / 180)
  ) +
  labs(
    title = "Near-bottom Data Availability",
    subtitle = "Taiwan Strait | Oct 2021 - Oct 2022",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  theme_bw()













