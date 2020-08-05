library(shiny)
library(tidyverse)
library(rnoaa)
library(plotly)
library(lubridate)
library(zoo)

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
                          limit = limit)}

# monitors <- c("USC00457933", "USW00024157", "USW00024114", "ASN00041023",
#               "ASN00009998", "ASN00066078", "ASN00003069", "ASN00090162",
#               "ASN00040126", "ASN00058161")
# obs <- meteo_pull_monitors(monitors)
# obs_covr <- meteo_coverage(obs)
# 
# if (interactive()) {
#     library("ggplot2")
#     autoplot(obs_covr)
# }

get_weather <- function(zipcode, start_date = as.Date('1961-01-01'), end_date = Sys.Date(), radius = 200, limit = 50){
    if(!is.Date(start_date)){start_date <- as.Date(start_date)}
    if(!is.Date(end_date)){end_date <- as.Date(end_date)}
    latLon <- get_zipcode_info(zipcode) %>% select(id, latitude, longitude)
    print(latLon)
    stations <- get_weather_stations(latLon, year(start_date), year(end_date), radius, limit)[[1]]
    print(stations)
    current_date <- start_date
    weather_master <- list()
    num_years_data <- 1
    while(current_date < end_date){
        i <- 0
        tmax <- c()
        tmax$data <- c()
        tmin <- c()
        tmin$data <- c()
        print(str_c("Gathering weather data for ", year(current_date)))
        while((i < length(stations$id)) && (length(tmax$data) == 0 || length(tmin$data) == 0)){
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
                         token = "HnmvXmMXFNeHpkLROUmJndwOyDPXATFJ")
            glimpse(tmax$data)
            glimpse(tmin$data)}

        if(length(tmax$data) == 0 || length(tmin$data) == 0){
            print("No weather data found, please increase limit or radius")
            return(weather_master)}else{
            print("Weather data found, recording now.")}
        weather_master[[num_years_data]] <- list(year = year(current_date),
                                                 station_id = stations$id[i],
                                                 tmax = tmax, 
                                                 tmin = tmin)
        num_years_data <- num_years_data + 1
        current_date <- current_date %m+% years(1)}
    return(weather_master)
}

# count <- 1 
# master <- list()
# zipcodes <- c("02111", "10001", "27514", "33101", "22434")
# for(i in zipcodes){
#     print(i)
#     master[[count]] <- list(zipcode = i,
#                             data = get_weather(i))
#     
#     count <- count + 1 
# }



#Convert TMIN/TMAX to TAVG for standardized temperature methodology
yakima_weather <- read.csv("./dat/yakima_washington.csv") 
yakima_weather$TAVG <- rowMeans(yakima_weather[5:6], na.rm = FALSE)
yakima_weather$TAVG.F <- yakima_weather$TAVG
yakima_weather$TAVG <- NULL
yakima_weather$TAVG.K <- sapply(yakima_weather$TAVG.F, function(x){return(((x - 32) * 5/9) + 273.15)})
yakima_weather$unicells <- sapply(yakima_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[1, "b0"],
                          m = metabolics[1, "mass_max"],
                          E = metabolics[1, "E"],
                          temp = x))})
yakima_weather$plants <- sapply(yakima_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[2, "b0"],
                          m = metabolics[2, "mass_max"],
                          E = metabolics[2, "E"],
                          temp = x))})
yakima_weather$invertebrates <- sapply(yakima_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[3, "b0"],
                          m = metabolics[3, "mass_max"],
                          E = metabolics[3, "E"],
                          temp = x))})
yakima_weather$amphibians <- sapply(yakima_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[4, "b0"],
                          m = metabolics[4, "mass_max"],
                          E = metabolics[4, "E"],
                          temp = x))})
yakima_weather$reptiles <- sapply(yakima_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[5, "b0"],
                          m = metabolics[5, "mass_max"],
                          E = metabolics[5, "E"],
                          temp = x))})
yakima_weather$average <- sapply(yakima_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[6, "b0"],
                          m = metabolics[6, "mass_max"],
                          E = metabolics[6, "E"],
                          temp = x))})
yakima_weather$DATE <- as.Date(yakima_weather$DATE)

#Weather available as TAVG for Brazil
brazil_weather <- read.csv("./dat/campinas_brazil.csv")
brazil_weather$TAVG.F <- brazil_weather$TAVG
brazil_weather$TAVG <- NULL
brazil_weather$TAVG.K <- sapply(brazil_weather$TAVG.F, function(x){return(((x - 32) * 5/9) + 273.15)})
brazil_weather$unicells <- sapply(brazil_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[1, "b0"],
                          m = metabolics[1, "mass_max"],
                          E = metabolics[1, "E"],
                          temp = x))})
brazil_weather$plants <- sapply(brazil_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[2, "b0"],
                          m = metabolics[2, "mass_max"],
                          E = metabolics[2, "E"],
                          temp = x))})
brazil_weather$invertebrates <- sapply(brazil_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[3, "b0"],
                          m = metabolics[3, "mass_max"],
                          E = metabolics[3, "E"],
                          temp = x))})
brazil_weather$amphibians <- sapply(brazil_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[4, "b0"],
                          m = metabolics[4, "mass_max"],
                          E = metabolics[4, "E"],
                          temp = x))})
brazil_weather$reptiles <- sapply(brazil_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[5, "b0"],
                          m = metabolics[5, "mass_max"],
                          E = metabolics[5, "E"],
                          temp = x))})
brazil_weather$average <- sapply(brazil_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[6, "b0"],
                          m = metabolics[6, "mass_max"],
                          E = metabolics[6, "E"],
                          temp = x))})
brazil_weather$DATE <- as.Date(brazil_weather$DATE)

#Weather available as TAVG for Brazil
mexico_weather <- read.csv("./dat/tlaxcala_mexico.csv")
mexico_weather$TAVG.F <- mexico_weather$TAVG
mexico_weather$TAVG <- NULL
mexico_weather$TAVG.K <- sapply(mexico_weather$TAVG.F, function(x){return(((x - 32) * 5/9) + 273.15)})
mexico_weather$unicells <- sapply(mexico_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[1, "b0"],
                          m = metabolics[1, "mass_max"],
                          E = metabolics[1, "E"],
                          temp = x))})
mexico_weather$plants <- sapply(mexico_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[2, "b0"],
                          m = metabolics[2, "mass_max"],
                          E = metabolics[2, "E"],
                          temp = x))})
mexico_weather$invertebrates <- sapply(mexico_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[3, "b0"],
                          m = metabolics[3, "mass_max"],
                          E = metabolics[3, "E"],
                          temp = x))})
mexico_weather$amphibians <- sapply(mexico_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[4, "b0"],
                          m = metabolics[4, "mass_max"],
                          E = metabolics[4, "E"],
                          temp = x))})
mexico_weather$reptiles <- sapply(mexico_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[5, "b0"],
                          m = metabolics[5, "mass_max"],
                          E = metabolics[5, "E"],
                          temp = x))})
mexico_weather$average <- sapply(mexico_weather$TAVG.K, function(x){
    return(metabolic_rate(b0 = metabolics[6, "b0"],
                          m = metabolics[6, "mass_max"],
                          E = metabolics[6, "E"],
                          temp = x))})
mexico_weather$DATE <- as.Date(mexico_weather$DATE)

#means are 5 year intervals beginning at the start date, with the exception of the first entry which is the standard reference period mean (1961-1990)
brazil_climate <- data.frame(
    start_date = c(
        as.Date('1961-01-01'), 
        as.Date('1980-01-01'), 
        as.Date('1985-01-01'), 
        as.Date('1990-01-01'), 
        as.Date('1995-01-01'), 
        as.Date('2000-01-01'), 
        as.Date('2005-01-01'), 
        as.Date('2010-01-01'), 
        as.Date('2015-01-01')),
    end_date = c(
        as.Date('1990-01-01'), 
        as.Date('1985-01-01'), 
        as.Date('1990-01-01'), 
        as.Date('1995-01-01'), 
        as.Date('2000-01-01'), 
        as.Date('2005-01-01'), 
        as.Date('2010-01-01'), 
        as.Date('2015-01-01'), 
        as.Date('2020-01-01')),
    mean_tavg = c(
        mean(brazil_weather$TAVG.K[which(brazil_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(brazil_weather$TAVG.K[which((brazil_weather$DATE >= as.Date('1980-01-01')) & (brazil_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(brazil_weather$TAVG.K[which((brazil_weather$DATE >= as.Date('1985-01-01')) & (brazil_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(brazil_weather$TAVG.K[which((brazil_weather$DATE >= as.Date('1990-01-01')) & (brazil_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(brazil_weather$TAVG.K[which((brazil_weather$DATE >= as.Date('1995-01-01')) & (brazil_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(brazil_weather$TAVG.K[which((brazil_weather$DATE >= as.Date('2000-01-01')) & (brazil_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(brazil_weather$TAVG.K[which((brazil_weather$DATE >= as.Date('2005-01-01')) & (brazil_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(brazil_weather$TAVG.K[which((brazil_weather$DATE >= as.Date('2010-01-01')) & (brazil_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(brazil_weather$TAVG.K[which((brazil_weather$DATE >= as.Date('2015-01-01')) & (brazil_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_reptiles = c(
        mean(brazil_weather$reptiles[which(brazil_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(brazil_weather$reptiles[which((brazil_weather$DATE >= as.Date('1980-01-01')) & (brazil_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(brazil_weather$reptiles[which((brazil_weather$DATE >= as.Date('1985-01-01')) & (brazil_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(brazil_weather$reptiles[which((brazil_weather$DATE >= as.Date('1990-01-01')) & (brazil_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(brazil_weather$reptiles[which((brazil_weather$DATE >= as.Date('1995-01-01')) & (brazil_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(brazil_weather$reptiles[which((brazil_weather$DATE >= as.Date('2000-01-01')) & (brazil_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(brazil_weather$reptiles[which((brazil_weather$DATE >= as.Date('2005-01-01')) & (brazil_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(brazil_weather$reptiles[which((brazil_weather$DATE >= as.Date('2010-01-01')) & (brazil_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(brazil_weather$reptiles[which((brazil_weather$DATE >= as.Date('2015-01-01')) & (brazil_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_unicells = c(
        mean(brazil_weather$unicells[which(brazil_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(brazil_weather$unicells[which((brazil_weather$DATE >= as.Date('1980-01-01')) & (brazil_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(brazil_weather$unicells[which((brazil_weather$DATE >= as.Date('1985-01-01')) & (brazil_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(brazil_weather$unicells[which((brazil_weather$DATE >= as.Date('1990-01-01')) & (brazil_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(brazil_weather$unicells[which((brazil_weather$DATE >= as.Date('1995-01-01')) & (brazil_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(brazil_weather$unicells[which((brazil_weather$DATE >= as.Date('2000-01-01')) & (brazil_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(brazil_weather$unicells[which((brazil_weather$DATE >= as.Date('2005-01-01')) & (brazil_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(brazil_weather$unicells[which((brazil_weather$DATE >= as.Date('2010-01-01')) & (brazil_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(brazil_weather$unicells[which((brazil_weather$DATE >= as.Date('2015-01-01')) & (brazil_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_invertebrates = c(
        mean(brazil_weather$invertebrates[which(brazil_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(brazil_weather$invertebrates[which((brazil_weather$DATE >= as.Date('1980-01-01')) & (brazil_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(brazil_weather$invertebrates[which((brazil_weather$DATE >= as.Date('1985-01-01')) & (brazil_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(brazil_weather$invertebrates[which((brazil_weather$DATE >= as.Date('1990-01-01')) & (brazil_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(brazil_weather$invertebrates[which((brazil_weather$DATE >= as.Date('1995-01-01')) & (brazil_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(brazil_weather$invertebrates[which((brazil_weather$DATE >= as.Date('2000-01-01')) & (brazil_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(brazil_weather$invertebrates[which((brazil_weather$DATE >= as.Date('2005-01-01')) & (brazil_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(brazil_weather$invertebrates[which((brazil_weather$DATE >= as.Date('2010-01-01')) & (brazil_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(brazil_weather$invertebrates[which((brazil_weather$DATE >= as.Date('2015-01-01')) & (brazil_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_amphibians = c(
        mean(brazil_weather$amphibians[which(brazil_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(brazil_weather$amphibians[which((brazil_weather$DATE >= as.Date('1980-01-01')) & (brazil_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(brazil_weather$amphibians[which((brazil_weather$DATE >= as.Date('1985-01-01')) & (brazil_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(brazil_weather$amphibians[which((brazil_weather$DATE >= as.Date('1990-01-01')) & (brazil_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(brazil_weather$amphibians[which((brazil_weather$DATE >= as.Date('1995-01-01')) & (brazil_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(brazil_weather$amphibians[which((brazil_weather$DATE >= as.Date('2000-01-01')) & (brazil_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(brazil_weather$amphibians[which((brazil_weather$DATE >= as.Date('2005-01-01')) & (brazil_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(brazil_weather$amphibians[which((brazil_weather$DATE >= as.Date('2010-01-01')) & (brazil_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(brazil_weather$amphibians[which((brazil_weather$DATE >= as.Date('2015-01-01')) & (brazil_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_plants = c(
        mean(brazil_weather$plants[which(brazil_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(brazil_weather$plants[which((brazil_weather$DATE >= as.Date('1980-01-01')) & (brazil_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(brazil_weather$plants[which((brazil_weather$DATE >= as.Date('1985-01-01')) & (brazil_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(brazil_weather$plants[which((brazil_weather$DATE >= as.Date('1990-01-01')) & (brazil_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(brazil_weather$plants[which((brazil_weather$DATE >= as.Date('1995-01-01')) & (brazil_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(brazil_weather$plants[which((brazil_weather$DATE >= as.Date('2000-01-01')) & (brazil_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(brazil_weather$plants[which((brazil_weather$DATE >= as.Date('2005-01-01')) & (brazil_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(brazil_weather$plants[which((brazil_weather$DATE >= as.Date('2010-01-01')) & (brazil_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(brazil_weather$plants[which((brazil_weather$DATE >= as.Date('2015-01-01')) & (brazil_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_average = c(
        mean(brazil_weather$average[which(brazil_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(brazil_weather$average[which((brazil_weather$DATE >= as.Date('1980-01-01')) & (brazil_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(brazil_weather$average[which((brazil_weather$DATE >= as.Date('1985-01-01')) & (brazil_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(brazil_weather$average[which((brazil_weather$DATE >= as.Date('1990-01-01')) & (brazil_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(brazil_weather$average[which((brazil_weather$DATE >= as.Date('1995-01-01')) & (brazil_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(brazil_weather$average[which((brazil_weather$DATE >= as.Date('2000-01-01')) & (brazil_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(brazil_weather$average[which((brazil_weather$DATE >= as.Date('2005-01-01')) & (brazil_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(brazil_weather$average[which((brazil_weather$DATE >= as.Date('2010-01-01')) & (brazil_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(brazil_weather$average[which((brazil_weather$DATE >= as.Date('2015-01-01')) & (brazil_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE))
)

#means are 5 year intervals beginning at the start date, with the exception of the first entry which is the standard reference period mean (1961-1990)
yakima_climate <- data.frame(
    start_date = c(
        as.Date('1961-01-01'), 
        as.Date('1980-01-01'), 
        as.Date('1985-01-01'), 
        as.Date('1990-01-01'), 
        as.Date('1995-01-01'), 
        as.Date('2000-01-01'), 
        as.Date('2005-01-01'), 
        as.Date('2010-01-01'), 
        as.Date('2015-01-01')),
    end_date = c(
        as.Date('1990-01-01'), 
        as.Date('1985-01-01'), 
        as.Date('1990-01-01'), 
        as.Date('1995-01-01'), 
        as.Date('2000-01-01'), 
        as.Date('2005-01-01'), 
        as.Date('2010-01-01'), 
        as.Date('2015-01-01'), 
        as.Date('2020-01-01')),
    mean_tavg = c(
        mean(yakima_weather$TAVG.K[which(yakima_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(yakima_weather$TAVG.K[which((yakima_weather$DATE >= as.Date('1980-01-01')) & (yakima_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(yakima_weather$TAVG.K[which((yakima_weather$DATE >= as.Date('1985-01-01')) & (yakima_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(yakima_weather$TAVG.K[which((yakima_weather$DATE >= as.Date('1990-01-01')) & (yakima_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(yakima_weather$TAVG.K[which((yakima_weather$DATE >= as.Date('1995-01-01')) & (yakima_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(yakima_weather$TAVG.K[which((yakima_weather$DATE >= as.Date('2000-01-01')) & (yakima_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(yakima_weather$TAVG.K[which((yakima_weather$DATE >= as.Date('2005-01-01')) & (yakima_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(yakima_weather$TAVG.K[which((yakima_weather$DATE >= as.Date('2010-01-01')) & (yakima_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(yakima_weather$TAVG.K[which((yakima_weather$DATE >= as.Date('2015-01-01')) & (yakima_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_reptiles = c(
        mean(yakima_weather$reptiles[which(yakima_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(yakima_weather$reptiles[which((yakima_weather$DATE >= as.Date('1980-01-01')) & (yakima_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(yakima_weather$reptiles[which((yakima_weather$DATE >= as.Date('1985-01-01')) & (yakima_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(yakima_weather$reptiles[which((yakima_weather$DATE >= as.Date('1990-01-01')) & (yakima_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(yakima_weather$reptiles[which((yakima_weather$DATE >= as.Date('1995-01-01')) & (yakima_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(yakima_weather$reptiles[which((yakima_weather$DATE >= as.Date('2000-01-01')) & (yakima_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(yakima_weather$reptiles[which((yakima_weather$DATE >= as.Date('2005-01-01')) & (yakima_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(yakima_weather$reptiles[which((yakima_weather$DATE >= as.Date('2010-01-01')) & (yakima_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(yakima_weather$reptiles[which((yakima_weather$DATE >= as.Date('2015-01-01')) & (yakima_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_unicells = c(
        mean(yakima_weather$unicells[which(yakima_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(yakima_weather$unicells[which((yakima_weather$DATE >= as.Date('1980-01-01')) & (yakima_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(yakima_weather$unicells[which((yakima_weather$DATE >= as.Date('1985-01-01')) & (yakima_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(yakima_weather$unicells[which((yakima_weather$DATE >= as.Date('1990-01-01')) & (yakima_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(yakima_weather$unicells[which((yakima_weather$DATE >= as.Date('1995-01-01')) & (yakima_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(yakima_weather$unicells[which((yakima_weather$DATE >= as.Date('2000-01-01')) & (yakima_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(yakima_weather$unicells[which((yakima_weather$DATE >= as.Date('2005-01-01')) & (yakima_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(yakima_weather$unicells[which((yakima_weather$DATE >= as.Date('2010-01-01')) & (yakima_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(yakima_weather$unicells[which((yakima_weather$DATE >= as.Date('2015-01-01')) & (yakima_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_invertebrates = c(
        mean(yakima_weather$invertebrates[which(yakima_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(yakima_weather$invertebrates[which((yakima_weather$DATE >= as.Date('1980-01-01')) & (yakima_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(yakima_weather$invertebrates[which((yakima_weather$DATE >= as.Date('1985-01-01')) & (yakima_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(yakima_weather$invertebrates[which((yakima_weather$DATE >= as.Date('1990-01-01')) & (yakima_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(yakima_weather$invertebrates[which((yakima_weather$DATE >= as.Date('1995-01-01')) & (yakima_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(yakima_weather$invertebrates[which((yakima_weather$DATE >= as.Date('2000-01-01')) & (yakima_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(yakima_weather$invertebrates[which((yakima_weather$DATE >= as.Date('2005-01-01')) & (yakima_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(yakima_weather$invertebrates[which((yakima_weather$DATE >= as.Date('2010-01-01')) & (yakima_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(yakima_weather$invertebrates[which((yakima_weather$DATE >= as.Date('2015-01-01')) & (yakima_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_amphibians = c(
        mean(yakima_weather$amphibians[which(yakima_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(yakima_weather$amphibians[which((yakima_weather$DATE >= as.Date('1980-01-01')) & (yakima_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(yakima_weather$amphibians[which((yakima_weather$DATE >= as.Date('1985-01-01')) & (yakima_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(yakima_weather$amphibians[which((yakima_weather$DATE >= as.Date('1990-01-01')) & (yakima_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(yakima_weather$amphibians[which((yakima_weather$DATE >= as.Date('1995-01-01')) & (yakima_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(yakima_weather$amphibians[which((yakima_weather$DATE >= as.Date('2000-01-01')) & (yakima_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(yakima_weather$amphibians[which((yakima_weather$DATE >= as.Date('2005-01-01')) & (yakima_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(yakima_weather$amphibians[which((yakima_weather$DATE >= as.Date('2010-01-01')) & (yakima_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(yakima_weather$amphibians[which((yakima_weather$DATE >= as.Date('2015-01-01')) & (yakima_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_plants = c(
        mean(yakima_weather$plants[which(yakima_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(yakima_weather$plants[which((yakima_weather$DATE >= as.Date('1980-01-01')) & (yakima_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(yakima_weather$plants[which((yakima_weather$DATE >= as.Date('1985-01-01')) & (yakima_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(yakima_weather$plants[which((yakima_weather$DATE >= as.Date('1990-01-01')) & (yakima_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(yakima_weather$plants[which((yakima_weather$DATE >= as.Date('1995-01-01')) & (yakima_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(yakima_weather$plants[which((yakima_weather$DATE >= as.Date('2000-01-01')) & (yakima_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(yakima_weather$plants[which((yakima_weather$DATE >= as.Date('2005-01-01')) & (yakima_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(yakima_weather$plants[which((yakima_weather$DATE >= as.Date('2010-01-01')) & (yakima_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(yakima_weather$plants[which((yakima_weather$DATE >= as.Date('2015-01-01')) & (yakima_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_average = c(
        mean(yakima_weather$average[which(yakima_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(yakima_weather$average[which((yakima_weather$DATE >= as.Date('1980-01-01')) & (yakima_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(yakima_weather$average[which((yakima_weather$DATE >= as.Date('1985-01-01')) & (yakima_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(yakima_weather$average[which((yakima_weather$DATE >= as.Date('1990-01-01')) & (yakima_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(yakima_weather$average[which((yakima_weather$DATE >= as.Date('1995-01-01')) & (yakima_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(yakima_weather$average[which((yakima_weather$DATE >= as.Date('2000-01-01')) & (yakima_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(yakima_weather$average[which((yakima_weather$DATE >= as.Date('2005-01-01')) & (yakima_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(yakima_weather$average[which((yakima_weather$DATE >= as.Date('2010-01-01')) & (yakima_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(yakima_weather$average[which((yakima_weather$DATE >= as.Date('2015-01-01')) & (yakima_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE))
)

mexico_climate <- data.frame(
    start_date = c(
        as.Date('1961-01-01'), 
        as.Date('1980-01-01'), 
        as.Date('1985-01-01'), 
        as.Date('1990-01-01'), 
        as.Date('1995-01-01'), 
        as.Date('2000-01-01'), 
        as.Date('2005-01-01'), 
        as.Date('2010-01-01'), 
        as.Date('2015-01-01')),
    end_date = c(
        as.Date('1990-01-01'), 
        as.Date('1985-01-01'), 
        as.Date('1990-01-01'), 
        as.Date('1995-01-01'), 
        as.Date('2000-01-01'), 
        as.Date('2005-01-01'), 
        as.Date('2010-01-01'), 
        as.Date('2015-01-01'), 
        as.Date('2020-01-01')),
    mean_tavg = c(
        mean(mexico_weather$TAVG.K[which(mexico_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(mexico_weather$TAVG.K[which((mexico_weather$DATE >= as.Date('1980-01-01')) & (mexico_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(mexico_weather$TAVG.K[which((mexico_weather$DATE >= as.Date('1985-01-01')) & (mexico_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(mexico_weather$TAVG.K[which((mexico_weather$DATE >= as.Date('1990-01-01')) & (mexico_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(mexico_weather$TAVG.K[which((mexico_weather$DATE >= as.Date('1995-01-01')) & (mexico_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(mexico_weather$TAVG.K[which((mexico_weather$DATE >= as.Date('2000-01-01')) & (mexico_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(mexico_weather$TAVG.K[which((mexico_weather$DATE >= as.Date('2005-01-01')) & (mexico_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(mexico_weather$TAVG.K[which((mexico_weather$DATE >= as.Date('2010-01-01')) & (mexico_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(mexico_weather$TAVG.K[which((mexico_weather$DATE >= as.Date('2015-01-01')) & (mexico_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_reptiles = c(
        mean(mexico_weather$reptiles[which(mexico_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(mexico_weather$reptiles[which((mexico_weather$DATE >= as.Date('1980-01-01')) & (mexico_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(mexico_weather$reptiles[which((mexico_weather$DATE >= as.Date('1985-01-01')) & (mexico_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(mexico_weather$reptiles[which((mexico_weather$DATE >= as.Date('1990-01-01')) & (mexico_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(mexico_weather$reptiles[which((mexico_weather$DATE >= as.Date('1995-01-01')) & (mexico_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(mexico_weather$reptiles[which((mexico_weather$DATE >= as.Date('2000-01-01')) & (mexico_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(mexico_weather$reptiles[which((mexico_weather$DATE >= as.Date('2005-01-01')) & (mexico_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(mexico_weather$reptiles[which((mexico_weather$DATE >= as.Date('2010-01-01')) & (mexico_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(mexico_weather$reptiles[which((mexico_weather$DATE >= as.Date('2015-01-01')) & (mexico_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_unicells = c(
        mean(mexico_weather$unicells[which(mexico_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(mexico_weather$unicells[which((mexico_weather$DATE >= as.Date('1980-01-01')) & (mexico_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(mexico_weather$unicells[which((mexico_weather$DATE >= as.Date('1985-01-01')) & (mexico_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(mexico_weather$unicells[which((mexico_weather$DATE >= as.Date('1990-01-01')) & (mexico_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(mexico_weather$unicells[which((mexico_weather$DATE >= as.Date('1995-01-01')) & (mexico_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(mexico_weather$unicells[which((mexico_weather$DATE >= as.Date('2000-01-01')) & (mexico_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(mexico_weather$unicells[which((mexico_weather$DATE >= as.Date('2005-01-01')) & (mexico_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(mexico_weather$unicells[which((mexico_weather$DATE >= as.Date('2010-01-01')) & (mexico_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(mexico_weather$unicells[which((mexico_weather$DATE >= as.Date('2015-01-01')) & (mexico_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_invertebrates = c(
        mean(mexico_weather$invertebrates[which(mexico_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(mexico_weather$invertebrates[which((mexico_weather$DATE >= as.Date('1980-01-01')) & (mexico_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(mexico_weather$invertebrates[which((mexico_weather$DATE >= as.Date('1985-01-01')) & (mexico_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(mexico_weather$invertebrates[which((mexico_weather$DATE >= as.Date('1990-01-01')) & (mexico_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(mexico_weather$invertebrates[which((mexico_weather$DATE >= as.Date('1995-01-01')) & (mexico_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(mexico_weather$invertebrates[which((mexico_weather$DATE >= as.Date('2000-01-01')) & (mexico_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(mexico_weather$invertebrates[which((mexico_weather$DATE >= as.Date('2005-01-01')) & (mexico_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(mexico_weather$invertebrates[which((mexico_weather$DATE >= as.Date('2010-01-01')) & (mexico_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(mexico_weather$invertebrates[which((mexico_weather$DATE >= as.Date('2015-01-01')) & (mexico_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_amphibians = c(
        mean(mexico_weather$amphibians[which(mexico_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(mexico_weather$amphibians[which((mexico_weather$DATE >= as.Date('1980-01-01')) & (mexico_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(mexico_weather$amphibians[which((mexico_weather$DATE >= as.Date('1985-01-01')) & (mexico_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(mexico_weather$amphibians[which((mexico_weather$DATE >= as.Date('1990-01-01')) & (mexico_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(mexico_weather$amphibians[which((mexico_weather$DATE >= as.Date('1995-01-01')) & (mexico_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(mexico_weather$amphibians[which((mexico_weather$DATE >= as.Date('2000-01-01')) & (mexico_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(mexico_weather$amphibians[which((mexico_weather$DATE >= as.Date('2005-01-01')) & (mexico_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(mexico_weather$amphibians[which((mexico_weather$DATE >= as.Date('2010-01-01')) & (mexico_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(mexico_weather$amphibians[which((mexico_weather$DATE >= as.Date('2015-01-01')) & (mexico_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_plants = c(
        mean(mexico_weather$plants[which(mexico_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(mexico_weather$plants[which((mexico_weather$DATE >= as.Date('1980-01-01')) & (mexico_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(mexico_weather$plants[which((mexico_weather$DATE >= as.Date('1985-01-01')) & (mexico_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(mexico_weather$plants[which((mexico_weather$DATE >= as.Date('1990-01-01')) & (mexico_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(mexico_weather$plants[which((mexico_weather$DATE >= as.Date('1995-01-01')) & (mexico_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(mexico_weather$plants[which((mexico_weather$DATE >= as.Date('2000-01-01')) & (mexico_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(mexico_weather$plants[which((mexico_weather$DATE >= as.Date('2005-01-01')) & (mexico_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(mexico_weather$plants[which((mexico_weather$DATE >= as.Date('2010-01-01')) & (mexico_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(mexico_weather$plants[which((mexico_weather$DATE >= as.Date('2015-01-01')) & (mexico_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
    mean_average = c(
        mean(mexico_weather$average[which(mexico_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(mexico_weather$average[which((mexico_weather$DATE >= as.Date('1980-01-01')) & (mexico_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(mexico_weather$average[which((mexico_weather$DATE >= as.Date('1985-01-01')) & (mexico_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(mexico_weather$average[which((mexico_weather$DATE >= as.Date('1990-01-01')) & (mexico_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(mexico_weather$average[which((mexico_weather$DATE >= as.Date('1995-01-01')) & (mexico_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(mexico_weather$average[which((mexico_weather$DATE >= as.Date('2000-01-01')) & (mexico_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(mexico_weather$average[which((mexico_weather$DATE >= as.Date('2005-01-01')) & (mexico_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(mexico_weather$average[which((mexico_weather$DATE >= as.Date('2010-01-01')) & (mexico_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(mexico_weather$average[which((mexico_weather$DATE >= as.Date('2015-01-01')) & (mexico_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE))
)

mexico_climate <- mexico_climate %>% 
    mutate(median_date = start_date + (as.numeric(difftime(end_date, start_date)) / 2))
brazil_climate <- brazil_climate %>% 
    mutate(median_date = start_date + (as.numeric(difftime(end_date, start_date)) / 2))
yakima_climate <- yakima_climate %>% 
    mutate(median_date = start_date + (as.numeric(difftime(end_date, start_date)) / 2))



plot_taxon <- function(taxon = "average"){
    column_id <- str_c("mean_", taxon)
    #climate (air temperature) plot
    #FIX START DATES (Point locations)
    climate_plot <- plot_ly(brazil_climate[2:length(brazil_climate[[1]]),], x = ~end_date) %>%
        #add_lines(y = ~TAVG.K) %>% 
        add_trace(data = yakima_climate[2:length(yakima_climate[[1]]),], y = ~mean_tavg, mode = 'lines+markers', name = "Washington, USA", color = 'green') %>%
        add_trace(data = brazil_climate[2:length(brazil_climate[[1]]),], y = ~mean_tavg, mode = 'lines+markers', name = "Campinas, BR", color = 'red') %>% 
        add_trace(data = mexico_climate[2:length(mexico_climate[[1]]),], y = ~mean_tavg, mode = 'lines+markers', name = "Tlaxcala, MX", color = 'blue') %>% 
        add_segments(x = mexico_climate[2, "start_date"], 
                     xend = mexico_climate[9, "end_date"], 
                     y = mexico_climate[1, "mean_tavg"], 
                     yend = mexico_climate[1, "mean_tavg"], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "green"),
                     showlegend = FALSE) %>% 
        add_segments(x = brazil_climate[2, "start_date"], 
                     xend = brazil_climate[9, "end_date"], 
                     y = brazil_climate[1, "mean_tavg"], 
                     yend = brazil_climate[1, "mean_tavg"], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "blue"),
                     showlegend = FALSE) %>% 
        add_segments(x = yakima_climate[2, "start_date"], 
                     xend = yakima_climate[9, "end_date"], 
                     y = yakima_climate[1, "mean_tavg"], 
                     yend = yakima_climate[1, "mean_tavg"], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "red")) %>% 
        layout(
            title = str_c("Climate means at 5 year intervals after reference period"),
            xaxis = list(
                title = "Date",
                rangeselector = list(
                    buttons = list(
                        list(
                            count = 1,
                            label = "1 y",
                            step = "year",
                            stepmode = "backward"),
                        list(
                            count = 6,
                            label = "6 y",
                            step = "year",
                            stepmode = "backward"),
                        list(
                            count = 20,
                            label = "20 y",
                            step = "year",
                            stepmode = "backward"))),
                
                rangeslider = list(type = "date")),
            yaxis = list(title = "Air Temperature \n(\u00B0K)"))
    
    
    metabolism_plot <- plot_ly(brazil_climate[2:length(brazil_climate[[1]]),], x = ~end_date) %>%
        #add_lines(y = ~TAVG.K) %>% 
        add_trace(data = yakima_climate[2:length(yakima_climate[[1]]),], y = yakima_climate[2:length(yakima_climate[[1]]),column_id], mode = 'lines+markers', name = "Washington, USA", color = 'green', showlegend = FALSE) %>%
        add_trace(data = brazil_climate[2:length(brazil_climate[[1]]),], y = brazil_climate[2:length(brazil_climate[[1]]),column_id], mode = 'lines+markers', name = "Campinas, BR", color = 'red', showlegend = FALSE) %>% 
        add_trace(data = mexico_climate[2:length(mexico_climate[[1]]),], y = mexico_climate[2:length(mexico_climate[[1]]),column_id], mode = 'lines+markers', name = "Tlaxcala, MX", color = 'blue', showlegend = FALSE) %>%
        add_segments(x = mexico_climate[2, "start_date"], 
                     xend = mexico_climate[9, "end_date"], 
                     y = mexico_climate[1, column_id], 
                     yend = mexico_climate[1, column_id], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "green"),
                     showlegend = FALSE) %>% 
        add_segments(x = brazil_climate[2, "start_date"], 
                     xend = brazil_climate[9, "end_date"], 
                     y = brazil_climate[1, column_id], 
                     yend = brazil_climate[1, column_id], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "blue"),
                     showlegend = FALSE) %>% 
        add_segments(x = yakima_climate[2, "start_date"], 
                     xend = yakima_climate[9, "end_date"], 
                     y = yakima_climate[1, column_id], 
                     yend = yakima_climate[1, column_id], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "red"),
                     showlegend = FALSE) %>% 
        layout(
            title = str_c("Metabolism and temperature means at 5 year intervals after reference period"),
            xaxis = list(
                title = "Date",
                rangeselector = list(
                    buttons = list(
                        list(
                            count = 1,
                            label = "1 y",
                            step = "year",
                            stepmode = "backward"),
                        list(
                            count = 6,
                            label = "6 y",
                            step = "year",
                            stepmode = "backward"),
                        list(
                            count = 20,
                            label = "20 y",
                            step = "year",
                            stepmode = "backward"))),
                
                rangeslider = list(type = "date")),
            yaxis = list(title = "Metabolic Rate \n(mWg<sup>-3/4</sup>)"))
    
    master_plots <- subplot(climate_plot, metabolism_plot, shareX = TRUE, nrows = 2, titleY = TRUE) 
    return(master_plots)}

# Define UI for application 
ui <- fluidPage(

    # Application title
    titlePanel("Climate Warming and Ectotherm Metabolism"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("taxon",
                        "Select a taxonomic group:",
                        choices = list("Unicells" = "unicells",
                                       "Plants" = "plants",
                                       "Invertebrates" = "invertebrates",
                                       "Amphibians" = "amphibians",
                                       "Reptiles" = "reptiles",
                                       "Average" = "average"),
                        selected = "average")
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotlyOutput("metaplot")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    #taxon <- reactive({input$taxon})
    #taxon <- "average"
    output$metaplot <- renderPlotly(plot_taxon(input$taxon))
}

# Run the application 
shinyApp(ui = ui, server = server)
