library(shiny)
library(tidyverse)
library(rnoaa)
library(plotly)
library(lubridate)
library(zoo)
library(shinyWidgets)
library(shinycssloaders)
library(shinytoastr)
library(mathjaxr)
# library(mapdeck)
# library(sp)
library(leaflet)
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


#Query and return relavant information about a zipcode
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
#Read
get_taxon_metabolics <- function(taxon, mass, weather_path = "./dat/yakima_washington.csv"){
    #Convert TMIN/TMAX to TAVG for standardized temperature methodology
    weather <- read.csv(weather_path) 
    weather$TAVG <- rowMeans(weather[5:6], na.rm = FALSE)
    weather$TAVG.F <- weather$TAVG
    weather$TAVG <- NULL
    weather$TAVG.K <- sapply(weather$TAVG.F, function(x){return(((x - 32) * 5/9) + 273.15)})
    weather$DATE <- as.Date(weather$DATE)
    # weather$unicells <- sapply(weather$TAVG.K, function(x){
    #     return(metabolic_rate(b0 = metabolics[1, "b0"],
    #                           m = (metabolics[1, "mass_min"] + metabolics[1, "mass_max"])/2,
    #                           E = metabolics[1, "E"],
    #                           temp = x))})
    weather[taxon] <- sapply(weather$TAVG.K, function(x){
        return(metabolic_rate(b0 = metabolics[1, "b0"],
                              m = mass,
                              E = metabolics[1, "E"],
                              temp = x))})
        
    return(weather)
}

get_taxon_averages <- function(taxon, taxon_metabolics){
    weather <- taxon_metabolics
    column_id <- str_c("mean_", taxon)
    climate <- data.frame(
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
            mean(weather$TAVG.K[which(weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
            mean(weather$TAVG.K[which((weather$DATE >= as.Date('1980-01-01')) & (weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
            mean(weather$TAVG.K[which((weather$DATE >= as.Date('1985-01-01')) & (weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
            mean(weather$TAVG.K[which((weather$DATE >= as.Date('1990-01-01')) & (weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
            mean(weather$TAVG.K[which((weather$DATE >= as.Date('1995-01-01')) & (weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
            mean(weather$TAVG.K[which((weather$DATE >= as.Date('2000-01-01')) & (weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
            mean(weather$TAVG.K[which((weather$DATE >= as.Date('2005-01-01')) & (weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
            mean(weather$TAVG.K[which((weather$DATE >= as.Date('2010-01-01')) & (weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
            mean(weather$TAVG.K[which((weather$DATE >= as.Date('2015-01-01')) & (weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)))
    climate[column_id] <- c(
        mean(weather[[taxon]][which(weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
        mean(weather[[taxon]][which((weather$DATE >= as.Date('1980-01-01')) & (weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
        mean(weather[[taxon]][which((weather$DATE >= as.Date('1985-01-01')) & (weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
        mean(weather[[taxon]][which((weather$DATE >= as.Date('1990-01-01')) & (weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
        mean(weather[[taxon]][which((weather$DATE >= as.Date('1995-01-01')) & (weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
        mean(weather[[taxon]][which((weather$DATE >= as.Date('2000-01-01')) & (weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
        mean(weather[[taxon]][which((weather$DATE >= as.Date('2005-01-01')) & (weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
        mean(weather[[taxon]][which((weather$DATE >= as.Date('2010-01-01')) & (weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
        mean(weather[[taxon]][which((weather$DATE >= as.Date('2015-01-01')) & (weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE))
    
    return(climate)
}


plot_taxon <- function(taxon = "average", mass){
    if(is.null(taxon)){taxon <- "average"}
    column_id <- str_c("mean_", taxon)
    mass <- as.integer(mass)
    available_mass <- ((metabolics[which(metabolics$taxon == taxon), "mass_min"] + metabolics[which(metabolics$taxon == taxon), "mass_max"])/2)
    available_mass <- as.integer(available_mass)
    if(mass == available_mass){
        toastr_info("Rendering default plots")
        mexico_climate <- read_rds(str_c("./dat/climate-analysis/", "tlaxcala_mexico", "_means.csv")) 
        mexico_weather <- read_rds(str_c("./dat/historical-weather/", "tlaxcala_mexico", ".csv")) 
        
        brazil_climate <- read_rds(str_c("./dat/climate-analysis/", "campinas_brazil", "_means.csv")) 
        brazil_weather <- read_rds(str_c("./dat/historical-weather/", "campinas_brazil", ".csv")) 
        
        yakima_climate <- read_rds(str_c("./dat/climate-analysis/", "yakima_washington", "_means.csv")) 
        yakima_weather <- read_rds(str_c("./dat/historical-weather/", "yakima_washington", ".csv")) 
        
        columbia_climate <- read_rds(str_c("./dat/climate-analysis/", "cali_alfonso_bonill_columbia", "_means.csv")) 
        columbia_weather <- read_rds(str_c("./dat/historical-weather/", "cali_alfonso_bonill_columbia", ".csv")) 
        
        greenland_climate <- read_rds(str_c("./dat/climate-analysis/", "danmarkshavn_greenland", "_means.csv")) 
        greenland_weather <- read_rds(str_c("./dat/historical-weather/", "danmarkshavn_greenland", ".csv")) 
        
    }else{
        toastr_info("Calculating metabolic parameters")
        mexico_weather <- get_taxon_metabolics(taxon = taxon, mass = mass, weather_path = str_c("./dat/", "tlaxcala_mexico", ".csv")) 
        mexico_climate <- get_taxon_averages(taxon = taxon, taxon_metabolics = mexico_weather)
        
        brazil_weather <- get_taxon_metabolics(taxon = taxon, mass = mass, weather_path = str_c("./dat/", "campinas_brazil", ".csv")) 
        brazil_climate <- get_taxon_averages(taxon = taxon, taxon_metabolics = brazil_weather)
        
        yakima_weather <- get_taxon_metabolics(taxon = taxon, mass = mass, weather_path = str_c("./dat/", "yakima_washington", ".csv")) 
        yakima_climate <- get_taxon_averages(taxon = taxon, taxon_metabolics = yakima_weather)
        
        columbia_weather <- get_taxon_metabolics(taxon = taxon, mass = mass, weather_path = str_c("./dat/", "cali_alfonso_bonill_columbia", ".csv")) 
        columbia_climate <- get_taxon_averages(taxon = taxon, taxon_metabolics = columbia_weather)
        
        greenland_weather <- get_taxon_metabolics(taxon = taxon, mass = mass, weather_path = str_c("./dat/", "danmarkshavn_greenland", ".csv")) 
        greenland_climate <- get_taxon_averages(taxon = taxon, taxon_metabolics = greenland_weather)}
    
    #climate (air temperature) plot
    #FIX START DATES (Point locations)
    climate_plot <- plot_ly(brazil_climate[2:length(brazil_climate[[1]]),], x = ~end_date) %>%
        #add_lines(y = ~TAVG.K) %>% 
        add_trace(type = "scatter", data = yakima_climate[2:length(yakima_climate[[1]]),], y = ~mean_tavg - yakima_climate[1,"mean_tavg"], mode = 'lines+markers', name = "Washington, USA", color = 'green') %>%
        add_trace(type = "scatter", data = brazil_climate[2:length(brazil_climate[[1]]),], y = ~mean_tavg - brazil_climate[1,"mean_tavg"], mode = 'lines+markers', name = "Campinas, BR", color = 'red') %>% 
        add_trace(type = "scatter", data = mexico_climate[2:length(mexico_climate[[1]]),], y = ~mean_tavg - mexico_climate[1,"mean_tavg"], mode = 'lines+markers', name = "Tlaxcala, MX", color = 'blue') %>% 
        add_trace(type = "scatter", data = columbia_climate[2:length(columbia_climate[[1]]),], y = ~mean_tavg - columbia_climate[1,"mean_tavg"], mode = 'lines+markers', name = "C. Alf. Bl., CO", color = 'pink') %>% 
        add_trace(type = "scatter", data = greenland_climate[2:length(greenland_climate[[1]]),], y = ~mean_tavg - greenland_climate[1,"mean_tavg"], mode = 'lines+markers', name = "Danmarkshavn, GL", color = 'violet') %>% 
        add_segments(x = greenland_climate[2, "start_date"], 
                     xend = greenland_climate[9, "end_date"], 
                     y = 0,
                     yend = 0,
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "grey"),
                     showlegend = TRUE) %>% 
        # add_segments(x = greenland_climate[2, "start_date"], 
        #              xend = greenland_climate[9, "end_date"], 
        #              y = greenland_climate[1, "mean_tavg"], 
        #              yend = greenland_climate[1, "mean_tavg"], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        # add_segments(x = mexico_climate[2, "start_date"], 
        #              xend = mexico_climate[9, "end_date"], 
        #              y = mexico_climate[1, "mean_tavg"], 
        #              yend = mexico_climate[1, "mean_tavg"], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        # add_segments(x = brazil_climate[2, "start_date"], 
        #              xend = brazil_climate[9, "end_date"], 
        #              y = brazil_climate[1, "mean_tavg"], 
        #              yend = brazil_climate[1, "mean_tavg"], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        # add_segments(x = yakima_climate[2, "start_date"], 
        #              xend = yakima_climate[9, "end_date"], 
        #              y = yakima_climate[1, "mean_tavg"], 
        #              yend = yakima_climate[1, "mean_tavg"], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey")) %>% 
        # add_segments(x = columbia_climate[2, "start_date"], 
        #              xend = columbia_climate[9, "end_date"], 
        #              y = columbia_climate[1, "mean_tavg"], 
        #              yend = columbia_climate[1, "mean_tavg"], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        layout(
            title = str_c("Climate means at 5 year intervals"),
            xaxis = list(
                showspikes = TRUE,
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
            yaxis = list(title = "Air Temp \n Change", ticksuffix = "\u00B0K"))
    
    
    metabolism_plot <- plot_ly(brazil_climate[2:length(brazil_climate[[1]]),], x = ~end_date) %>%
        #add_lines(y = ~TAVG.K) %>% 
        add_trace(type = "scatter", data = yakima_climate[2:length(yakima_climate[[1]]),], y = yakima_climate[2:length(yakima_climate[[1]]),column_id] - yakima_climate[1,column_id], mode = 'lines+markers', name = "Washington, USA", color = 'green', showlegend = FALSE) %>%
        add_trace(type = "scatter", data = brazil_climate[2:length(brazil_climate[[1]]),], y = brazil_climate[2:length(brazil_climate[[1]]),column_id] - brazil_climate[1,column_id], mode = 'lines+markers', name = "Campinas, BR", color = 'red', showlegend = FALSE) %>% 
        add_trace(type = "scatter", data = mexico_climate[2:length(mexico_climate[[1]]),], y = mexico_climate[2:length(mexico_climate[[1]]),column_id] - mexico_climate[1,column_id], mode = 'lines+markers', name = "Tlaxcala, MX", color = 'blue', showlegend = FALSE) %>%
        add_trace(type = "scatter", data = columbia_climate[2:length(columbia_climate[[1]]),], y = columbia_climate[2:length(columbia_climate[[1]]),column_id] - columbia_climate[1,column_id], mode = 'lines+markers', name = "C. Alf. Bl., CO", color = 'pink', showlegend = FALSE) %>%
        add_trace(type = "scatter", data = greenland_climate[2:length(greenland_climate[[1]]),], y = greenland_climate[2:length(greenland_climate[[1]]),column_id] - greenland_climate[1,column_id], mode = 'lines+markers', name = "Danmarkshavn, GL", color = 'violet', showlegend = FALSE) %>%
        add_segments(x = greenland_climate[2, "start_date"], 
                     xend = greenland_climate[9, "end_date"], 
                     y = 0,
                     yend = 0,
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "grey"),
                     showlegend = FALSE) %>% 
        # add_segments(x = greenland_climate[2, "start_date"], 
        #              xend = greenland_climate[9, "end_date"], 
        #              y = greenland_climate[1, column_id], 
        #              yend = greenland_climate[1, column_id], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        # add_segments(x = mexico_climate[2, "start_date"], 
        #              xend = mexico_climate[9, "end_date"], 
        #              y = mexico_climate[1, column_id], 
        #              yend = mexico_climate[1, column_id], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        # add_segments(x = brazil_climate[2, "start_date"], 
        #              xend = brazil_climate[9, "end_date"], 
        #              y = brazil_climate[1, column_id], 
        #              yend = brazil_climate[1, column_id], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        # add_segments(x = yakima_climate[2, "start_date"], 
        #              xend = yakima_climate[9, "end_date"], 
        #              y = yakima_climate[1, column_id], 
        #              yend = yakima_climate[1, column_id], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        # add_segments(x = columbia_climate[2, "start_date"], 
        #              xend = columbia_climate[9, "end_date"], 
        #              y = columbia_climate[1, column_id], 
        #              yend = columbia_climate[1, column_id], 
        #              name = "Reference Period",
        #              line = list(dash = "dot",
        #                          color = "grey"),
        #              showlegend = FALSE) %>% 
        layout(
            title = str_c("Metabolism and temperature means at 5 year intervals"),
            xaxis = list(
                title = "Date",
                showspikes = TRUE,
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
            yaxis = list(title = "Metabolic \nRate Change \n(mWg<sup>-3/4</sup>)"))
    
    
    
    metabolism_percent_plot <- plot_ly(brazil_climate[2:length(brazil_climate[[1]]),], x = ~end_date) %>%
        #add_lines(y = ~TAVG.K) %>%
        add_trace(type = "scatter", data = yakima_climate[2:length(yakima_climate[[1]]),], y = ((yakima_climate[2:length(yakima_climate[[1]]),column_id] - yakima_climate[1, column_id]) / yakima_climate[1, column_id])*100, mode = 'lines+markers', name = "Washington, USA", color = 'green', showlegend = FALSE) %>%
        add_trace(type = "scatter", data = brazil_climate[2:length(brazil_climate[[1]]),], y = ((brazil_climate[2:length(brazil_climate[[1]]),column_id] - brazil_climate[1, column_id]) / brazil_climate[1, column_id])*100, mode = 'lines+markers', name = "Campinas, BR", color = 'red', showlegend = FALSE) %>%
        add_trace(type = "scatter", data = mexico_climate[2:length(mexico_climate[[1]]),], y = ((mexico_climate[2:length(mexico_climate[[1]]),column_id] - mexico_climate[1, column_id]) / mexico_climate[1, column_id])*100, mode = 'lines+markers', name = "Tlaxcala, MX", color = 'blue', showlegend = FALSE) %>%
        add_trace(type = "scatter", data = columbia_climate[2:length(columbia_climate[[1]]),], y = ((columbia_climate[2:length(columbia_climate[[1]]),column_id] - columbia_climate[1, column_id]) / columbia_climate[1, column_id])*100, mode = 'lines+markers', name = "C. Alf. Bl., CO", color = 'pink', showlegend = FALSE) %>%
        add_trace(type = "scatter", data = greenland_climate[2:length(greenland_climate[[1]]),], y = ((greenland_climate[2:length(greenland_climate[[1]]),column_id] - greenland_climate[1, column_id]) / greenland_climate[1, column_id])*100, mode = 'lines+markers', name = "Danmarkshavn, GL", color = 'violet', showlegend = FALSE) %>%
        add_segments(x = greenland_climate[2, "start_date"], 
                     xend = greenland_climate[9, "end_date"], 
                     y = greenland_climate[1, column_id], 
                     yend = greenland_climate[1, column_id], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "grey"),
                     showlegend = FALSE) %>% 
        add_segments(x = mexico_climate[2, "start_date"], 
                     xend = mexico_climate[9, "end_date"], 
                     y = mexico_climate[1, column_id], 
                     yend = mexico_climate[1, column_id], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "grey"),
                     showlegend = FALSE) %>% 
        add_segments(x = brazil_climate[2, "start_date"], 
                     xend = brazil_climate[9, "end_date"], 
                     y = brazil_climate[1, column_id], 
                     yend = brazil_climate[1, column_id], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "grey"),
                     showlegend = FALSE) %>% 
        add_segments(x = yakima_climate[2, "start_date"], 
                     xend = yakima_climate[9, "end_date"], 
                     y = yakima_climate[1, column_id], 
                     yend = yakima_climate[1, column_id], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "grey"),
                     showlegend = FALSE) %>% 
        add_segments(x = columbia_climate[2, "start_date"], 
                     xend = columbia_climate[9, "end_date"], 
                     y = columbia_climate[1, column_id], 
                     yend = columbia_climate[1, column_id], 
                     name = "Reference Period",
                     line = list(dash = "dot",
                                 color = "grey"),
                     showlegend = FALSE) %>% 
        layout(
            title = str_c("Metabolism and temperature means at 5 year intervals"),
            xaxis = list(
                title = "Date",
                showspikes = TRUE,
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
            yaxis = list(title = "Metabolic \nRate Change", ticksuffix = "%"))
    
    master_plots <- subplot(climate_plot, metabolism_plot, metabolism_percent_plot, shareX = TRUE, nrows = 3, titleY = TRUE) 
    toastr_success("Plots rendered")
    return(master_plots)}

# Define UI for application 
ui <- fluidPage(
    # Application title
    useToastr(),
    withMathJax(),
    tags$head(tags$link(rel="shortcut icon", href="./favicon.ico")),
    titlePanel("Climate Warming and Ectotherm Metabolism"),
    hr(),
    # Sidebar with a slider input for number of bins 
    verticalLayout(
        #includeMarkdown("./intro1.md"),
        includeMarkdown("./finintro.md"),
        leafletOutput("climate_map") %>% withSpinner(color = "#228B22"),
        helpText("The selected weather stations' (black markers) positions span multiple climate zones and represent a large latitudinal range from the equator to the north pole."),
        # includeMarkdown("./intro2.md"),
        includeMarkdown("./finintro2.md"),
        # p("$$B(m,T) = b_{0}m^{3/4}e^{-E/kT}$$"),
        # includeMarkdown("./intro3.md"),
        #hr(),
        dropdownButton(
            tags$h3("Plot Settings"),
            circle = TRUE, status = "info", icon = icon("gear"), width = "300px",
            tooltip = tooltipOptions(title = "Click to see plot settings!"),
            #bg_color = "#c4d5a7", width = 2L,
            #bg_color = "#dec4de", color = "#228B22",
            selectInput("taxon",
                        "Taxonomic group:",
                        choices = list("Unicells" = "unicells",
                                       "Plants" = "plants",
                                       "Invertebrates" = "invertebrates",
                                       "Amphibians" = "amphibians",
                                       "Reptiles" = "reptiles",
                                       "Average" = "average"),
                        selected = "average"),
            sliderInput(
                inputId = "mass",
                label = "Taxon mass",
                value = (metabolics$mass_max[6] + metabolics$mass_min[6])/2,
                min = metabolics$mass_min[6],
                max = metabolics$mass_max[6],
                post = "g",
                animate = FALSE#, 
                # lineCap = "round",
                # fgColor = "#428BCA",
                # inputColor = "#428BCA",
                # width = "150px",
                # immediate = FALSE
            )),
        # Show a plot of the generated distribution
        
        plotlyOutput("metaplot", height = "100%") %>% withSpinner(color = "#228B22"),
        helpText("Air Temperature means, metabolic rate means, 
                 and the percent metabolic rate change from baseline (1961-1990, 
                 standard reference period) are plotted here for each weather station. 
                 Grey dashed lines represent means of the standard reference period.
                 Click the settings icon above to change the taxon and mass inputs 
                 in this interactive graphic. ")
    ))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    #taxon <- reactive({input$taxon})
    
    observeEvent(input$taxon, {
        if(!is.null(input$taxon)){
            metabolic_inf <- metabolics[which(metabolics$taxon == input$taxon),]
            updateSliderInput(
                session = session,
                inputId = "mass",
                value = (metabolic_inf$mass_max + metabolic_inf$mass_min)/2,
                min = metabolic_inf$mass_min,
                max = metabolic_inf$mass_max,
                step = ceiling({(metabolic_inf$mass_max - metabolic_inf$mass_min)/10}),
            )
            output$metaplot <- renderPlotly(plot_taxon(taxon = input$taxon, mass = (metabolic_inf$mass_max + metabolic_inf$mass_min)/2))}
    })
    
    observeEvent(input$mass, {
        metabolic_inf <- metabolics[which(metabolics$taxon == input$taxon),]
        if(as.integer(input$mass) != as.integer((metabolic_inf$mass_max + metabolic_inf$mass_min)/2)){
        output$metaplot <- renderPlotly(plot_taxon(taxon = input$taxon, mass = input$mass))}
    })
    
    colors = c("red", "orange", "green", "blue")
    zones = c("Tropical zone", "Subtropical zone", "Temperate zone", "Polar/subpolar zone")
    # pal <- colorBin(colors, 
    #                 zones, 
    #                 bins = zones,
    #                 na.color = "transparent")
    #CartoDB.DarkMatter
    climate_map <- leaflet() %>%
        addProviderTiles(providers$CartoDB.Positron, options = providerTileOptions(noWrap = FALSE)) %>% 
        flyTo(0,0,1) %>% 
        addRectangles(180, 66.5, -180, 90, stroke = FALSE, color = "blue", fillOpacity = 0.4) %>% 
        addRectangles(180, -66.5, -180, -90, stroke = FALSE, color = "blue", fillOpacity = 0.4) %>% 
        addRectangles(180, 23.5, -180, 35, stroke = FALSE, color = "orange", fillOpacity = 0.4) %>% 
        addRectangles(180, -23.5, -180, -35, stroke = FALSE, color = "orange", fillOpacity = 0.4) %>% 
        addRectangles(180, 35, -180, 66.5, stroke = FALSE, color = "green", fillOpacity = 0.4) %>% 
        addRectangles(180, -35, -180, -66.5, stroke = FALSE, color = "green", fillOpacity = 0.4) %>% 
        addRectangles(180, -23.5, -180, 23.5, stroke = FALSE, color = "red", fillOpacity = 0.4) %>% 
        addCircleMarkers(lat = 46.6021, lng = -120.5059, color = "black", radius = 7, 
                         label = "Yakima, Washington, USA") %>% 
        addCircleMarkers(lat = -22.9329, lng = -47.0738, color = "black", radius = 7,
                         label = "Campinas, Brazil") %>% 
        addCircleMarkers(lat = 19.3182, lng = -98.2375, color = "black", radius = 7,
                         label = "Tlaxcala, Mexico") %>% 
        addCircleMarkers(lat = 3.5411, lng = -76.3846, color = "black", radius = 7,
                         label = "Cali Alfonso Bonill, Columbia") %>% 
        addCircleMarkers(lat = 76.7667, lng = -18.6667, color = "black", radius = 7,
                         label = "Danmarkshavn, Greenland") %>% 
        addLegend(position = "bottomright", 
                  colors = colors, 
                  labels = zones,
                  title = "Climate Zones") 
    
    output$climate_map <- renderLeaflet(climate_map) 
    
}

# Run the application 
shinyApp(ui = ui, server = server)
