getwd()
setwd("C:/Users/Mostafa/Desktop/Review project/BFAST/20")

# Load necessary libraries
library(sp)
library(raster)
library(bfast)
library(gdalUtilities)
library(sf)
library(foreach)
library(doParallel)
library(tictoc)  # Optional library for timing the operations

# Read your TIFF file containing the time series NBR data
tiff_file <- "NBR15072024.tif"
nbr_ts <- brick(tiff_file)
dates <- read.csv("Date.csv")

# Assign dates to time-series NBR images
image_dates <- as.Date(dates$date, format = "%Y-%m-%d")
names(nbr_ts) <- as.character(image_dates)

# Define the start year for monitoring
start_monitoring_year <- 2004

# Function to calculate the number of disturbances for a single pixel
calculate_disturbances <- function(pixel_values, image_dates, start_monitoring_year) {
  # Replace zero with NA
  pixel_values[pixel_values == 0] <- NA
  
  # Calculate Frequency
  time_diff <- diff(image_dates)
  avg_time_diff <- mean(time_diff)
  avg_time_diff_days <- as.numeric(avg_time_diff, units = "days")
  avg_time_diff_months <- avg_time_diff_days / 30
  frequency <- round(12 / avg_time_diff_months)
  
  # Convert NBR values to time series object
  nbr_val_ts <- ts(pixel_values, start = start(image_dates), frequency = frequency)
  
  # Apply BFAST algorithm
  bfast_result <- tryCatch({
    bfast(nbr_val_ts, max.iter = 10, h = 10)
  }, error = function(e) {
    return(NULL)
  })
  
  if (is.null(bfast_result) || is.null(bfast_result$output[[1]]$bp.Vt)) {
    return(0)
  }
  
  # Extracting breakpoint information
  breakpoints_info <- bfast_result$output[[1]]$bp.Vt
  
  # Check if breakpoints_info is a list and contains breakpoints
  if (!is.list(breakpoints_info) || is.null(breakpoints_info$breakpoints)) {
    return(0)
  }
  
  breakpoint_obs_numbers <- breakpoints_info$breakpoints
  
  # Return the number of breakpoints (disturbances)
  return(length(breakpoint_obs_numbers))
}

# Initialize a raster to store the number of disturbances
disturbances_raster <- raster(nbr_ts)

# Set up parallel processing
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Loop through each pixel and calculate the number of disturbances in parallel
tic("Parallel computation")
disturbances <- foreach(i = 1:ncell(nbr_ts), .combine = 'c', .packages = c('raster', 'bfast')) %dopar% {
  pixel_values <- as.numeric(nbr_ts[i])
  calculate_disturbances(pixel_values, image_dates, start_monitoring_year)
}
toc()

# Stop the cluster
stopCluster(cl)

# Assign the disturbances to the raster
disturbances_raster[] <- disturbances

# Frequency of trend
disturbance_count <- table(values(disturbances_raster))

# Function to retrieve coordinates
get_coordinates <- function(raster_layer, value) {
  idx <- which(values(raster_layer) == value)
  xy <- xyFromCell(raster_layer, idx)
  return(data.frame(x = xy[,1], y = xy[,2], value = value))
}

# Get coordinates for each disturbance count
unique_disturbance_values <- unique(values(disturbances_raster))
coords_list <- lapply(unique_disturbance_values, function(val) {
  get_coordinates(disturbances_raster, val)
})

# Convert UTM to Decimal
convert_utm_to_decimal <- function(coords, utm_crs) {
  coords_sp <- SpatialPoints(coords, proj4string = CRS(utm_crs))
  coords_latlon <- spTransform(coords_sp, CRS("+proj=longlat +datum=WGS84"))
  return(data.frame(longitude = coordinates(coords_latlon)[,1], latitude = coordinates(coords_latlon)[,2]))
}

# Get UTM CRS from the disturbances raster
utm_crs <- crs(disturbances_raster)

# Get coordinates for each disturbance count and convert them to decimal
unique_disturbance_values <- unique(values(disturbances_raster))
coords_list <- lapply(unique_disturbance_values, function(val) {
  coords <- get_coordinates(disturbances_raster, val)
  decimal_coords <- convert_utm_to_decimal(coords[, c('x', 'y')], utm_crs@projargs)
  decimal_coords$value <- val  # Add the disturbance value to the dataframe
  return(decimal_coords)
})

# Combine all coordinates into one dataframe
coordinates_df <- do.call(rbind, coords_list)
write.csv(coordinates_df, "disturbance_coordinates17072024.csv", row.names = FALSE)

# Set custom margins for the plot
par(mar = c(5, 5, 4, 2) + 0.1)
# Plot the number of disturbances with axis labels and legend title
plot(disturbances_raster, main = "Number of Detected Disturbances", col = rev(terrain.colors(13)),
     xlab = "X", ylab = "Y")
