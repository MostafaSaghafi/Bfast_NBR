# Load necessary libraries
library(shiny)
library(raster)
library(bfast)
library(sp)
library(RColorBrewer)
library(sf)

# Set maximum request size to 1024MB (you can adjust this value as needed)
options(shiny.maxRequestSize = 50*1024^20)

# Define UI
ui <- fluidPage(
  titlePanel("BFAST Algorithm on Raster Map"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("tiffFile", "Choose a TIFF file", accept = ".tif"),
      fileInput("dateFile", "Choose a Date CSV file", accept = ".csv"),
      numericInput("startYear", "Start Monitoring Year:", value = 2004, min = 1900, max = 2100),
      actionButton("process", "Process"),
      textInput("manualCoords", "Enter Coordinates (latitude, longitude):", value = ""),
      actionButton("processManual", "Process Manual Coordinates")
    ),
    
    mainPanel(
      plotOutput("rasterPlot", click = "rasterClick"),
      verbatimTextOutput("clickInfo"),
      plotOutput("bfastPlot")
    )
  )
)

# Define Server
server <- function(input, output) {
  rasterData <- reactiveValues(nbr_ts = NULL, dates = NULL, image_dates = NULL)
  
  # Observe TIFF file input
  observeEvent(input$tiffFile, {
    req(input$tiffFile)
    tiff_file <- input$tiffFile$datapath
    rasterData$nbr_ts <- brick(tiff_file)
  })
  
  # Observe Date CSV file input
  observeEvent(input$dateFile, {
    req(input$dateFile)
    dates <- read.csv(input$dateFile$datapath)
    rasterData$dates <- dates
    rasterData$image_dates <- as.Date(dates$date, format = "%Y-%m-%d")
    names(rasterData$nbr_ts) <- as.character(rasterData$image_dates)
  })
  
  # Render the raster plot
  output$rasterPlot <- renderPlot({
    req(rasterData$nbr_ts)
    
    if (nlayers(rasterData$nbr_ts) >= 3) {
      # Choose suitable bands for RGB visualization
      plotRGB(rasterData$nbr_ts, r = 3, g = 2, b = 1, stretch = "hist")
    } else {
      plot(rasterData$nbr_ts)
    }
  })
  
  # Print click information with coordinates in decimal degrees
  output$clickInfo <- renderPrint({
    req(input$rasterClick)
    coords <- as.numeric(c(input$rasterClick$x, input$rasterClick$y))
    
    # Convert the coordinates to decimal degrees
    raster_crs <- st_crs(rasterData$nbr_ts) # Get CRS from raster
    utm_coords <- st_as_sf(data.frame(x = coords[1], y = coords[2]), coords = c("x", "y"), crs = raster_crs)
    latlon_coords <- st_transform(utm_coords, crs = 4326) # Transform to WGS84
    
    print(st_coordinates(latlon_coords))
  })
  
  # Function to process coordinates
  processCoordinates <- function(coords, crs = 4326) {
    tryCatch({
      # Create sf object with the provided coordinates
      coords_sf <- st_as_sf(data.frame(x = coords[1], y = coords[2]), coords = c("x", "y"), crs = crs)
      
      # Transform coordinates to the raster's CRS
      raster_crs <- st_crs(rasterData$nbr_ts) # Get CRS from raster
      coords_transformed <- st_transform(coords_sf, crs = raster_crs)
      
      if (is.null(coords_transformed) || nrow(coords_transformed) == 0) {
        stop("Invalid coordinates after transformation")
      }
      
      stable_forest_pixel <- cellFromXY(rasterData$nbr_ts, st_coordinates(coords_transformed))
      if (is.na(stable_forest_pixel)) {
        stop("Invalid pixel coordinates.")
      }
      
      print(paste("Pixel Number:", stable_forest_pixel))
      
      nbr_val <- extract(rasterData$nbr_ts, stable_forest_pixel)
      nbr_val <- as.numeric(nbr_val)
      nbr_val[nbr_val == 0] <- NA
      
      if (length(nbr_val) == 0 || all(is.na(nbr_val))) {
        stop("No valid NBR values found.")
      }
      
      print("NBR Values:")
      print(nbr_val)
      
      time_diff <- diff(rasterData$image_dates)
      avg_time_diff <- mean(time_diff)
      avg_time_diff_days <- as.numeric(avg_time_diff, units = "days")
      avg_time_diff_months <- avg_time_diff_days / 30
      frequency <- round(12 / avg_time_diff_months)
      
      nbr_val_ts <- ts(nbr_val, start = c(input$startYear, 1), frequency = frequency)
      print("Time Series Data:")
      print(nbr_val_ts)
      
      bfast_result <- bfast(nbr_val_ts, max.iter = 10, h = 10)
      
      output$bfastPlot <- renderPlot({
        par(mfrow = c(1,1))
        plot(bfast_result, type = "components")
      })
      
      breakpoints_info <- bfast_result$output[[1]]$bp.Vt
      breakpoint_obs_numbers <- breakpoints_info$breakpoints
      breakpoint_dates <- rasterData$image_dates[breakpoint_obs_numbers]
      
      print("Breakpoint Observation Numbers:")
      print(breakpoint_obs_numbers)
      print("Breakpoint Dates:")
      print(breakpoint_dates)
    }, error = function(e) {
      print(paste("Error:", e$message))
    })
  }
  
  # Observe clicks on the raster plot
  observeEvent(input$rasterClick, {
    req(input$rasterClick)
    req(rasterData$nbr_ts)
    
    coords <- as.numeric(c(input$rasterClick$x, input$rasterClick$y))
    print(paste("Clicked coordinates:", coords[1], coords[2]))
    
    processCoordinates(coords, crs = st_crs(rasterData$nbr_ts))
  })
  
  # Observe manual coordinate input
  observeEvent(input$processManual, {
    req(input$manualCoords)
    req(rasterData$nbr_ts)
    
    coords <- as.numeric(unlist(strsplit(input$manualCoords, ",")))
    if (length(coords) != 2) {
      print("Invalid coordinates format. Please enter as 'latitude, longitude'.")
      return()
    }
    
    print(paste("Manual coordinates:", coords[1], coords[2]))
    
    processCoordinates(coords)
  })
}

# Run App
shinyApp(ui, server)