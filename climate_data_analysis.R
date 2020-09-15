library(tidyverse)
library(rnoaa)
library(lubridate)
library(zoo)
# library(mapdeck)
# library(sp)
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
  #return(b0*(m^(3/4))*exp(-E/(k*temp)))
  return(b0*exp(-E/(k*temp)))
}

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

#Convert TMIN/TMAX to TAVG for standardized temperature methodology
greenland_weather <- read.csv("./dat/danmarkshavn_greenland.csv") 
greenland_weather$TAVG <- rowMeans(greenland_weather[5:6], na.rm = FALSE)
greenland_weather$TAVG.F <- greenland_weather$TAVG
greenland_weather$TAVG <- NULL
greenland_weather$TAVG.K <- sapply(greenland_weather$TAVG.F, function(x){return(((x - 32) * 5/9) + 273.15)})
greenland_weather$unicells <- sapply(greenland_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[1, "b0"],
                        m = metabolics[1, "mass_max"],
                        E = metabolics[1, "E"],
                        temp = x))})
greenland_weather$plants <- sapply(greenland_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[2, "b0"],
                        m = metabolics[2, "mass_max"],
                        E = metabolics[2, "E"],
                        temp = x))})
greenland_weather$invertebrates <- sapply(greenland_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[3, "b0"],
                        m = metabolics[3, "mass_max"],
                        E = metabolics[3, "E"],
                        temp = x))})
greenland_weather$amphibians <- sapply(greenland_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[4, "b0"],
                        m = metabolics[4, "mass_max"],
                        E = metabolics[4, "E"],
                        temp = x))})
greenland_weather$reptiles <- sapply(greenland_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[5, "b0"],
                        m = metabolics[5, "mass_max"],
                        E = metabolics[5, "E"],
                        temp = x))})
greenland_weather$average <- sapply(greenland_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[6, "b0"],
                        m = metabolics[6, "mass_max"],
                        E = metabolics[6, "E"],
                        temp = x))})
greenland_weather$DATE <- as.Date(greenland_weather$DATE)

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

#Weather available as TAVG for Brazil
columbia_weather <- read.csv("./dat/cali_alfonso_bonill_columbia.csv")
columbia_weather$TAVG.F <- columbia_weather$TAVG
columbia_weather$TAVG <- NULL
columbia_weather$TAVG.K <- sapply(columbia_weather$TAVG.F, function(x){return(((x - 32) * 5/9) + 273.15)})
columbia_weather$unicells <- sapply(columbia_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[1, "b0"],
                        m = metabolics[1, "mass_max"],
                        E = metabolics[1, "E"],
                        temp = x))})
columbia_weather$plants <- sapply(columbia_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[2, "b0"],
                        m = metabolics[2, "mass_max"],
                        E = metabolics[2, "E"],
                        temp = x))})
columbia_weather$invertebrates <- sapply(columbia_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[3, "b0"],
                        m = metabolics[3, "mass_max"],
                        E = metabolics[3, "E"],
                        temp = x))})
columbia_weather$amphibians <- sapply(columbia_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[4, "b0"],
                        m = metabolics[4, "mass_max"],
                        E = metabolics[4, "E"],
                        temp = x))})
columbia_weather$reptiles <- sapply(columbia_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[5, "b0"],
                        m = metabolics[5, "mass_max"],
                        E = metabolics[5, "E"],
                        temp = x))})
columbia_weather$average <- sapply(columbia_weather$TAVG.K, function(x){
  return(metabolic_rate(b0 = metabolics[6, "b0"],
                        m = metabolics[6, "mass_max"],
                        E = metabolics[6, "E"],
                        temp = x))})
columbia_weather$
  columbia_weather$DATE <- as.Date(columbia_weather$DATE)

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

#means are 5 year intervals beginning at the start date, with the exception of the first entry which is the standard reference period mean (1961-1990)
greenland_climate <- data.frame(
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
    mean(greenland_weather$TAVG.K[which(greenland_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(greenland_weather$TAVG.K[which((greenland_weather$DATE >= as.Date('1980-01-01')) & (greenland_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(greenland_weather$TAVG.K[which((greenland_weather$DATE >= as.Date('1985-01-01')) & (greenland_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(greenland_weather$TAVG.K[which((greenland_weather$DATE >= as.Date('1990-01-01')) & (greenland_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(greenland_weather$TAVG.K[which((greenland_weather$DATE >= as.Date('1995-01-01')) & (greenland_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(greenland_weather$TAVG.K[which((greenland_weather$DATE >= as.Date('2000-01-01')) & (greenland_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(greenland_weather$TAVG.K[which((greenland_weather$DATE >= as.Date('2005-01-01')) & (greenland_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(greenland_weather$TAVG.K[which((greenland_weather$DATE >= as.Date('2010-01-01')) & (greenland_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(greenland_weather$TAVG.K[which((greenland_weather$DATE >= as.Date('2015-01-01')) & (greenland_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_reptiles = c(
    mean(greenland_weather$reptiles[which(greenland_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(greenland_weather$reptiles[which((greenland_weather$DATE >= as.Date('1980-01-01')) & (greenland_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(greenland_weather$reptiles[which((greenland_weather$DATE >= as.Date('1985-01-01')) & (greenland_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(greenland_weather$reptiles[which((greenland_weather$DATE >= as.Date('1990-01-01')) & (greenland_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(greenland_weather$reptiles[which((greenland_weather$DATE >= as.Date('1995-01-01')) & (greenland_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(greenland_weather$reptiles[which((greenland_weather$DATE >= as.Date('2000-01-01')) & (greenland_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(greenland_weather$reptiles[which((greenland_weather$DATE >= as.Date('2005-01-01')) & (greenland_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(greenland_weather$reptiles[which((greenland_weather$DATE >= as.Date('2010-01-01')) & (greenland_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(greenland_weather$reptiles[which((greenland_weather$DATE >= as.Date('2015-01-01')) & (greenland_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_unicells = c(
    mean(greenland_weather$unicells[which(greenland_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(greenland_weather$unicells[which((greenland_weather$DATE >= as.Date('1980-01-01')) & (greenland_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(greenland_weather$unicells[which((greenland_weather$DATE >= as.Date('1985-01-01')) & (greenland_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(greenland_weather$unicells[which((greenland_weather$DATE >= as.Date('1990-01-01')) & (greenland_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(greenland_weather$unicells[which((greenland_weather$DATE >= as.Date('1995-01-01')) & (greenland_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(greenland_weather$unicells[which((greenland_weather$DATE >= as.Date('2000-01-01')) & (greenland_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(greenland_weather$unicells[which((greenland_weather$DATE >= as.Date('2005-01-01')) & (greenland_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(greenland_weather$unicells[which((greenland_weather$DATE >= as.Date('2010-01-01')) & (greenland_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(greenland_weather$unicells[which((greenland_weather$DATE >= as.Date('2015-01-01')) & (greenland_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_invertebrates = c(
    mean(greenland_weather$invertebrates[which(greenland_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(greenland_weather$invertebrates[which((greenland_weather$DATE >= as.Date('1980-01-01')) & (greenland_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(greenland_weather$invertebrates[which((greenland_weather$DATE >= as.Date('1985-01-01')) & (greenland_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(greenland_weather$invertebrates[which((greenland_weather$DATE >= as.Date('1990-01-01')) & (greenland_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(greenland_weather$invertebrates[which((greenland_weather$DATE >= as.Date('1995-01-01')) & (greenland_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(greenland_weather$invertebrates[which((greenland_weather$DATE >= as.Date('2000-01-01')) & (greenland_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(greenland_weather$invertebrates[which((greenland_weather$DATE >= as.Date('2005-01-01')) & (greenland_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(greenland_weather$invertebrates[which((greenland_weather$DATE >= as.Date('2010-01-01')) & (greenland_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(greenland_weather$invertebrates[which((greenland_weather$DATE >= as.Date('2015-01-01')) & (greenland_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_amphibians = c(
    mean(greenland_weather$amphibians[which(greenland_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(greenland_weather$amphibians[which((greenland_weather$DATE >= as.Date('1980-01-01')) & (greenland_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(greenland_weather$amphibians[which((greenland_weather$DATE >= as.Date('1985-01-01')) & (greenland_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(greenland_weather$amphibians[which((greenland_weather$DATE >= as.Date('1990-01-01')) & (greenland_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(greenland_weather$amphibians[which((greenland_weather$DATE >= as.Date('1995-01-01')) & (greenland_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(greenland_weather$amphibians[which((greenland_weather$DATE >= as.Date('2000-01-01')) & (greenland_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(greenland_weather$amphibians[which((greenland_weather$DATE >= as.Date('2005-01-01')) & (greenland_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(greenland_weather$amphibians[which((greenland_weather$DATE >= as.Date('2010-01-01')) & (greenland_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(greenland_weather$amphibians[which((greenland_weather$DATE >= as.Date('2015-01-01')) & (greenland_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_plants = c(
    mean(greenland_weather$plants[which(greenland_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(greenland_weather$plants[which((greenland_weather$DATE >= as.Date('1980-01-01')) & (greenland_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(greenland_weather$plants[which((greenland_weather$DATE >= as.Date('1985-01-01')) & (greenland_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(greenland_weather$plants[which((greenland_weather$DATE >= as.Date('1990-01-01')) & (greenland_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(greenland_weather$plants[which((greenland_weather$DATE >= as.Date('1995-01-01')) & (greenland_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(greenland_weather$plants[which((greenland_weather$DATE >= as.Date('2000-01-01')) & (greenland_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(greenland_weather$plants[which((greenland_weather$DATE >= as.Date('2005-01-01')) & (greenland_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(greenland_weather$plants[which((greenland_weather$DATE >= as.Date('2010-01-01')) & (greenland_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(greenland_weather$plants[which((greenland_weather$DATE >= as.Date('2015-01-01')) & (greenland_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_average = c(
    mean(greenland_weather$average[which(greenland_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(greenland_weather$average[which((greenland_weather$DATE >= as.Date('1980-01-01')) & (greenland_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(greenland_weather$average[which((greenland_weather$DATE >= as.Date('1985-01-01')) & (greenland_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(greenland_weather$average[which((greenland_weather$DATE >= as.Date('1990-01-01')) & (greenland_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(greenland_weather$average[which((greenland_weather$DATE >= as.Date('1995-01-01')) & (greenland_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(greenland_weather$average[which((greenland_weather$DATE >= as.Date('2000-01-01')) & (greenland_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(greenland_weather$average[which((greenland_weather$DATE >= as.Date('2005-01-01')) & (greenland_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(greenland_weather$average[which((greenland_weather$DATE >= as.Date('2010-01-01')) & (greenland_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(greenland_weather$average[which((greenland_weather$DATE >= as.Date('2015-01-01')) & (greenland_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE))
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

columbia_climate <- data.frame(
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
    mean(columbia_weather$TAVG.K[which(columbia_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(columbia_weather$TAVG.K[which((columbia_weather$DATE >= as.Date('1980-01-01')) & (columbia_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(columbia_weather$TAVG.K[which((columbia_weather$DATE >= as.Date('1985-01-01')) & (columbia_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(columbia_weather$TAVG.K[which((columbia_weather$DATE >= as.Date('1990-01-01')) & (columbia_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(columbia_weather$TAVG.K[which((columbia_weather$DATE >= as.Date('1995-01-01')) & (columbia_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(columbia_weather$TAVG.K[which((columbia_weather$DATE >= as.Date('2000-01-01')) & (columbia_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(columbia_weather$TAVG.K[which((columbia_weather$DATE >= as.Date('2005-01-01')) & (columbia_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(columbia_weather$TAVG.K[which((columbia_weather$DATE >= as.Date('2010-01-01')) & (columbia_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(columbia_weather$TAVG.K[which((columbia_weather$DATE >= as.Date('2015-01-01')) & (columbia_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_reptiles = c(
    mean(columbia_weather$reptiles[which(columbia_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(columbia_weather$reptiles[which((columbia_weather$DATE >= as.Date('1980-01-01')) & (columbia_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(columbia_weather$reptiles[which((columbia_weather$DATE >= as.Date('1985-01-01')) & (columbia_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(columbia_weather$reptiles[which((columbia_weather$DATE >= as.Date('1990-01-01')) & (columbia_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(columbia_weather$reptiles[which((columbia_weather$DATE >= as.Date('1995-01-01')) & (columbia_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(columbia_weather$reptiles[which((columbia_weather$DATE >= as.Date('2000-01-01')) & (columbia_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(columbia_weather$reptiles[which((columbia_weather$DATE >= as.Date('2005-01-01')) & (columbia_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(columbia_weather$reptiles[which((columbia_weather$DATE >= as.Date('2010-01-01')) & (columbia_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(columbia_weather$reptiles[which((columbia_weather$DATE >= as.Date('2015-01-01')) & (columbia_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_unicells = c(
    mean(columbia_weather$unicells[which(columbia_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(columbia_weather$unicells[which((columbia_weather$DATE >= as.Date('1980-01-01')) & (columbia_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(columbia_weather$unicells[which((columbia_weather$DATE >= as.Date('1985-01-01')) & (columbia_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(columbia_weather$unicells[which((columbia_weather$DATE >= as.Date('1990-01-01')) & (columbia_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(columbia_weather$unicells[which((columbia_weather$DATE >= as.Date('1995-01-01')) & (columbia_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(columbia_weather$unicells[which((columbia_weather$DATE >= as.Date('2000-01-01')) & (columbia_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(columbia_weather$unicells[which((columbia_weather$DATE >= as.Date('2005-01-01')) & (columbia_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(columbia_weather$unicells[which((columbia_weather$DATE >= as.Date('2010-01-01')) & (columbia_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(columbia_weather$unicells[which((columbia_weather$DATE >= as.Date('2015-01-01')) & (columbia_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_invertebrates = c(
    mean(columbia_weather$invertebrates[which(columbia_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(columbia_weather$invertebrates[which((columbia_weather$DATE >= as.Date('1980-01-01')) & (columbia_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(columbia_weather$invertebrates[which((columbia_weather$DATE >= as.Date('1985-01-01')) & (columbia_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(columbia_weather$invertebrates[which((columbia_weather$DATE >= as.Date('1990-01-01')) & (columbia_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(columbia_weather$invertebrates[which((columbia_weather$DATE >= as.Date('1995-01-01')) & (columbia_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(columbia_weather$invertebrates[which((columbia_weather$DATE >= as.Date('2000-01-01')) & (columbia_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(columbia_weather$invertebrates[which((columbia_weather$DATE >= as.Date('2005-01-01')) & (columbia_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(columbia_weather$invertebrates[which((columbia_weather$DATE >= as.Date('2010-01-01')) & (columbia_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(columbia_weather$invertebrates[which((columbia_weather$DATE >= as.Date('2015-01-01')) & (columbia_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_amphibians = c(
    mean(columbia_weather$amphibians[which(columbia_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(columbia_weather$amphibians[which((columbia_weather$DATE >= as.Date('1980-01-01')) & (columbia_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(columbia_weather$amphibians[which((columbia_weather$DATE >= as.Date('1985-01-01')) & (columbia_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(columbia_weather$amphibians[which((columbia_weather$DATE >= as.Date('1990-01-01')) & (columbia_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(columbia_weather$amphibians[which((columbia_weather$DATE >= as.Date('1995-01-01')) & (columbia_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(columbia_weather$amphibians[which((columbia_weather$DATE >= as.Date('2000-01-01')) & (columbia_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(columbia_weather$amphibians[which((columbia_weather$DATE >= as.Date('2005-01-01')) & (columbia_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(columbia_weather$amphibians[which((columbia_weather$DATE >= as.Date('2010-01-01')) & (columbia_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(columbia_weather$amphibians[which((columbia_weather$DATE >= as.Date('2015-01-01')) & (columbia_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_plants = c(
    mean(columbia_weather$plants[which(columbia_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(columbia_weather$plants[which((columbia_weather$DATE >= as.Date('1980-01-01')) & (columbia_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(columbia_weather$plants[which((columbia_weather$DATE >= as.Date('1985-01-01')) & (columbia_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(columbia_weather$plants[which((columbia_weather$DATE >= as.Date('1990-01-01')) & (columbia_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(columbia_weather$plants[which((columbia_weather$DATE >= as.Date('1995-01-01')) & (columbia_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(columbia_weather$plants[which((columbia_weather$DATE >= as.Date('2000-01-01')) & (columbia_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(columbia_weather$plants[which((columbia_weather$DATE >= as.Date('2005-01-01')) & (columbia_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(columbia_weather$plants[which((columbia_weather$DATE >= as.Date('2010-01-01')) & (columbia_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(columbia_weather$plants[which((columbia_weather$DATE >= as.Date('2015-01-01')) & (columbia_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)),
  mean_average = c(
    mean(columbia_weather$average[which(columbia_weather$DATE < as.Date('1990-01-01'))], na.rm = TRUE),
    mean(columbia_weather$average[which((columbia_weather$DATE >= as.Date('1980-01-01')) & (columbia_weather$DATE < as.Date('1985-01-01')))], na.rm = TRUE),
    mean(columbia_weather$average[which((columbia_weather$DATE >= as.Date('1985-01-01')) & (columbia_weather$DATE < as.Date('1990-01-01')))], na.rm = TRUE),
    mean(columbia_weather$average[which((columbia_weather$DATE >= as.Date('1990-01-01')) & (columbia_weather$DATE < as.Date('1995-01-01')))], na.rm = TRUE),
    mean(columbia_weather$average[which((columbia_weather$DATE >= as.Date('1995-01-01')) & (columbia_weather$DATE < as.Date('2000-01-01')))], na.rm = TRUE),
    mean(columbia_weather$average[which((columbia_weather$DATE >= as.Date('2000-01-01')) & (columbia_weather$DATE < as.Date('2005-01-01')))], na.rm = TRUE),
    mean(columbia_weather$average[which((columbia_weather$DATE >= as.Date('2005-01-01')) & (columbia_weather$DATE < as.Date('2010-01-01')))], na.rm = TRUE),
    mean(columbia_weather$average[which((columbia_weather$DATE >= as.Date('2010-01-01')) & (columbia_weather$DATE < as.Date('2015-01-01')))], na.rm = TRUE),
    mean(columbia_weather$average[which((columbia_weather$DATE >= as.Date('2015-01-01')) & (columbia_weather$DATE < as.Date('2020-01-01')))], na.rm = TRUE)))


mexico_climate <- mexico_climate %>% 
  mutate(median_date = start_date + (as.numeric(difftime(end_date, start_date)) / 2))
brazil_climate <- brazil_climate %>% 
  mutate(median_date = start_date + (as.numeric(difftime(end_date, start_date)) / 2))
yakima_climate <- yakima_climate %>% 
  mutate(median_date = start_date + (as.numeric(difftime(end_date, start_date)) / 2))
columbia_climate <- columbia_climate %>% 
  mutate(median_date = start_date + (as.numeric(difftime(end_date, start_date)) / 2))


write_rds(mexico_climate, str_c("./dat/climate-analysis/", "tlaxcala_mexico", "_means.csv")) 
write_rds(mexico_weather, str_c("./dat/historical-weather/", "tlaxcala_mexico", ".csv")) 

write_rds(brazil_climate, str_c("./dat/climate-analysis/", "campinas_brazil", "_means.csv")) 
write_rds(brazil_weather, str_c("./dat/historical-weather/", "campinas_brazil", ".csv")) 

write_rds(yakima_climate, str_c("./dat/climate-analysis/", "yakima_washington", "_means.csv")) 
write_rds(yakima_weather, str_c("./dat/historical-weather/", "yakima_washington", ".csv")) 

write_rds(columbia_climate, str_c("./dat/climate-analysis/", "cali_alfonso_bonill_columbia", "_means.csv")) 
write_rds(columbia_weather, str_c("./dat/historical-weather/", "cali_alfonso_bonill_columbia", ".csv")) 

write_rds(greenland_climate, str_c("./dat/climate-analysis/", "danmarkshavn_greenland", "_means.csv")) 
write_rds(greenland_weather, str_c("./dat/historical-weather/", "danmarkshavn_greenland", ".csv")) 

