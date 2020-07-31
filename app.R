#

library(shiny)
library(tidyverse)
library(rnoaa)
library(plotly)
library(lubridate)

#Data source: https://www.nature.com/articles/nature09407
#Mass is is grams 
metabolics <- data.frame("taxon" = c("unicells", "plants", "invertebrates", "amphibians", "reptiles", "average"),
                         "slope" = c(-8.79, -7.61, -9.15, -5.76, -8.78, -8.02),
                         "b0" = c(25.8, 21.37, 27.62, 16.68, 26.85, 23.66),
                         "E" = c(0.76, 0.66, 0.79, 0.5, 0.76, 0.69),
                         "mass_min" = c(0.000000000001, 0.1, 0.0001, 0.1, 10, 0.000000000001),
                         "mass_max" = c(0.1, 100, 100, 100, 10000, 10000))
# b0 --- an empirically derived and taxon specific normalization constant
# m ---- body mass 
# E ---- average activation energy for biochemical reactions of metabolism 
# temp - body temperature (K) (assumed to be same as air temperature for thermoconforming and exposed ectotherms)
metabolic_rate <- function(b0, m, E, temp){
    #k - Boltzmann constant (eV/K, relates temperature to energy) 
    k = 0.000086173
    return(b0*(m^(3/4))*exp(-E/(k*temp)))
}

get_zipcode_info <- function(zipcode){
    zipcodes <- read_delim("./dat/us-zip-code-latitude-and-longitude.csv", delim = ";")
    zipLocation <- zipcodes %>% 
        filter(Zip == zipcode) %>% 
        select(City, State, geopoint)
    print(str_c("Zip location: ", zipLocation$geopoint))
    if(!identical(zipLocation$geopoint, character(0))){
        
        latLon <- data.frame("id" = zipcode,
                       "city" = zipLocation$City,
                       "state" = zipLocation$State,
                       "latitude" = unlist(strsplit(zipLocation$geopoint, ","))[1],
                       "longitude" = unlist(strsplit(zipLocation$geopoint, ","))[2])
    }
    return(latLon)}

#---------------Only run the following if you want to update ghcnd-stations.txt------------
#This will query NOAA for the most up to date info on weather stations in the GHCND network
#The output is saved to a file called ghcnd-stations-current.csv
#in the specified "path/to/directory", or in the current working directory if no directory is specified

updateGHCNDStations <- function(directory = getwd()){
    print("Getting ghcnd-stations.txt from NOAA...")
    stationsDailyRaw <- read.fwf(url("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-stations.txt"),
                                 widths = c(11, 9, 11, 7, 2, 31, 5, 10),
                                 header = FALSE, strip.white = TRUE, comment.char = "",
                                 stringsAsFactors = FALSE)
    print("Getting ghcnd-inventory.txt from NOAA...")
    inventoryDailyRaw <- read.fwf(url("https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/ghcnd-inventory.txt"),
                                  widths = c(11, 9, 10, 5, 5, 5),
                                  header = FALSE, strip.white = TRUE, comment.char = "",
                                  stringsAsFactors = FALSE)
    stationColNames <- c("id","latitude", "longitude", "elevation",
                         "state", "name", "gsn_flag", "wmo_id")
    inventoryColNames <- c("id","latitude", "longitude",
                           "element", "first_year", "last_year")
    
    ghcndStationsDaily <- stats::setNames(stationsDailyRaw, stationColNames)
    ghcndInventoryDaily <- stats::setNames(inventoryDailyRaw, inventoryColNames)
    
    ghcndStationsDailyComplete <- merge(ghcndStationsDaily, ghcndInventoryDaily[, -c(2, 3)], by = "id")
    
    sturdyGHCNDStations <- tibble::as_tibble(ghcndStationsDailyComplete[stats::complete.cases(ghcndStationsDailyComplete), ])
    
    saveRDS(sturdyGHCNDStations, file = str_c(directory, "/ghcnd-stations-current.csv"))
    
    return(sturdyGHCNDStations)}

#--------End Update of GHCND Stations----------

#read in current local copy of ghcnd stations 
localGHCNDStations <- readRDS(file = "./dat/ghcnd-stations-current.csv")


# lat_lon_df must be a list with columns "id", "latitude", and "longitude"
# Example usage: 
# latLon <- get_zipcode_info(99203) %>% select(id, latitude, longitude)
# stations <- get_weather_stations(latLon, 1961, 2020, 500, 5)[[1]]
get_weather_stations <- function(lat_lon_df, year_min, year_max, radius, limit){
    meteo_nearby_stations(lat_lon_df = lat_lon_df,
                          station_data = localGHCNDStations,
                          var = c("TMAX", "TMIN"),
                          year_min = year_min,
                          year_max = year_max,
                          radius = radius,
                          limit = limit)
}

monitors <- c("ASN00095063", "ASN00024025", "ASN00040112", "ASN00041023",
              "ASN00009998", "ASN00066078", "ASN00003069", "ASN00090162",
              "ASN00040126", "ASN00058161")
obs <- meteo_pull_monitors(head(stations[1], 10)$id)
obs_covr <- meteo_coverage(obs)

if (interactive()) {
    library("ggplot2")
    autoplot(obs_covr)
}

get_weather <- function(zipcode, radius = 250, start_date = as.Date('1961-01-01'), end_date = Sys.Date()){
    if(!is.Date(start_date)){start_date <- as.Date(start_date)}
    if(!is.Date(end_date)){end_date <- as.Date(end_date)}
    latLon <- get_zipcode_info(zipcode) %>% select(id, latitude, longitude)
    print(latLon)
    stations <- get_weather_stations(latLon, year(start_date), year(end_date), radius, 30)[[1]]
    print(stations)
    current_date <- start_date
    weather_master <- list()
    num_years_data <- 1
    while(current_date < end_date){
        i <- 0
        tmax <- c()
        tmax$data[[1]] <- c()
        tmin <- c()
        tmin$data[[1]] <- c()
        print(str_c("Gathering weather data for ", year(current_date), " from ", stations$id[i]))
        while((length(tmax$data[[1]]) == 0 || length(tmin$data[[1]] == 0)) && i < length(stations$id)){
            if((start_date %m+% years(1)) < end_date)
                {temp_end_date <- start_date %m+% years(1)}
            else{temp_end_date <- end_date}
            i <- i + 1 
            print("Fetching TMAX")
            tmax <- ncdc(datasetid = "GHCND", 
                            stationid = paste0("GHCND:", stations$id[i]), 
                            datatypeid = "TMAX", 
                            startdate = current_date, 
                            enddate = temp_end_date, 
                            token = "HnmvXmMXFNeHpkLROUmJndwOyDPXATFJ")
            print("Fetching TMIN")
            tmin <- ncdc(datasetid = "GHCND", 
                         stationid = paste0("GHCND:", stations$id[i]), 
                         datatypeid = "TMAX", 
                         startdate = current_date, 
                         enddate = temp_end_date, 
                         token = "HnmvXmMXFNeHpkLROUmJndwOyDPXATFJ")}

        weather_master[[num_years_data]] <- list(year = year(current_date),
                                                 station_id = stations$id[i],
                                                 tmax = tmax, 
                                                 tmin = tmin)
        num_years_data <- num_years_data + 1
        current_date <- current_date %m+% years(1)}
    return(weather_master)
}


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Climate Warming and Ectotherm Metabolism"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("bins",
                        "Number of bins:",
                        min = 1,
                        max = 50,
                        value = 30)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
