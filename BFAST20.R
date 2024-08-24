getwd()
setwd("C:/Users/Mostafa/Desktop/BFAST/20")

# Load necessary libraries
library(raster)
library(bfast)
library(gdalUtilities)
# Read your TIFF file containing the time series NBR data
tiff_file <- "nbr20.tif"
nbr_ts <- brick(tiff_file)
dates <- read.csv("Date.csv")

# Assign dates to time-series NBR images
image_dates <- as.Date(dates$date, format = "%Y-%m-%d")
names(nbr_ts) <- as.character(image_dates)

# Define the start year for monitoring
start_monitoring_year <- 2004

# Get pixel values on specified coordinates (X=Long, Y=Lat)
lat_lon_coords <- SpatialPoints(cbind(46.4733, 35.5888), proj4string = CRS("+proj=longlat +datum=WGS84"))
utm_coords <- spTransform(lat_lon_coords, CRS("+proj=utm +zone=38 +datum=WGS84 +units=m +no_defs"))
stable_forest_pixel <- cellFromXY(nbr_ts, coordinates(utm_coords))

# Pull out NBR time series for a single pixel 
nbr_val <- as.numeric(nbr_ts[stable_forest_pixel])

# Replace zero with NA
nbr_val[nbr_val == 0] <- NA
plot(nbr_val)

##Calculate Frequency
# Calculate the time differences between consecutive images
time_diff <- diff(image_dates)
# Calculate the average time difference
avg_time_diff <- mean(time_diff)
# Convert the average time difference to days
avg_time_diff_days <- as.numeric(avg_time_diff, units = "days")
# Convert the average time difference to months (assuming 30 days per month)
avg_time_diff_months <- avg_time_diff_days / 30
# Determine the appropriate frequency based on the average time difference
frequency <- round(12 / avg_time_diff_months)
# Display the calculated frequency
print(frequency)

# Convert NBR values to time series object
nbr_val_ts <- ts(nbr_val, start = start(image_dates), frequency=frequency)

# Apply BFAST algorithm
bfast_result <- bfast(nbr_val_ts, max.iter = 10 , h=10)  # Adjust parameters as needed)
plot(bfast_result, type = "components")
plot(bfast_result, type = "all")
plot(bfast_result, type = "data")
plot(bfast_result, type = "seasonal")
plot(bfast_result, type = "trend")
plot(bfast_result, type = "noise")

##Extract BP date
# Extracting breakpoint information
breakpoints_info <- bfast_result$output[[1]]$bp.Vt
# Extracting observation number of breakpoints
breakpoint_obs_numbers <- breakpoints_info$breakpoints
breakpoint_dates <- image_dates[breakpoint_obs_numbers]
# Printing breakpoint observation numbers and corresponding dates
print(breakpoint_obs_numbers)
print(breakpoint_dates)
