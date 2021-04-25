# 24 November 2020
# Snapshot Wisconsin data cannot be published to protect volunteer privacy
# specifically, lat-longs associated with cameras cannot be published
# So here's code to simulate time-of-day detection data,
# fit it with a nonparametric kernel density model, 
# and calculate the proportion of activity within night, dawn, day, dusk

library(here)
library(sf)
library(tidyverse)
library(lubridate)
library(suncalc)
library(activity)

setwd(here::here("data"))

# shapefile of Wisconsin...for generating random points
wi <- sf::st_read("WI_state_outline.shp") %>% 
  sf::st_transform(crs = 4326)

# how many random points to simulate
npoints <- 10
# how many detections per site to simulate
ndets <- 10

# generate random points within Wisconsin
random_points <- sf::st_sample(wi, size = npoints) %>% 
  sf::st_as_sf(.) %>%
  dplyr::mutate(XCOORD = sf::st_coordinates(.)[,1], 
                YCOORD = sf::st_coordinates(.)[,2]) %>% 
  sf::st_drop_geometry(.)

# sanity check - do we get random points showing up within Wisconsin?
ggplot() + 
  geom_sf(data = wi, aes(geometry = geometry)) + 
  geom_point(data = random_points, aes(x = XCOORD, y = YCOORD))

# start date
start <- lubridate::ymd("2017-12-15")
# end date
end <- lubridate::ymd("2017-12-20")

dates <- seq(from = start,
             to = end,
             by = 1)

# create a dataframe with the unique points & dates of survey 
points_dates <- as_tibble(random_points[rep(1:nrow(random_points),
                                    times = length(dates)),]) %>% 
  arrange(XCOORD, YCOORD) %>% 
  add_column(date = rep(dates, npoints))

# create a dataframe with the points and simulated detections
random_detections <- as_tibble(random_points[rep(1:nrow(random_points),
                                                 times = ndets),]) %>% 
  arrange(XCOORD, YCOORD) %>% 
  # here, we're simulating detection date-times, 10 per site
  add_column(detection = as_datetime(start)
             + runif(ndets*npoints, 0, as.numeric(difftime(end,
                                                           start,
                                                           unit = "sec")))) %>% 
  mutate(date = as_date(detection))

# get times of relevant times of day for defining periods for each site
suntimes <- suncalc::getSunlightTimes(data = dplyr::select(points_dates, 
                                                           date, 
                                                           lat = YCOORD,
                                                           lon = XCOORD),
                                      keep = c("nightEnd",
                                               "goldenHourEnd",
                                               "goldenHour", 
                                               "night")) %>% 
  mutate_at(vars(nightEnd:night), function(x) with_tz(x, "America/Chicago")) %>% 
  rename(YCOORD = lat,
         XCOORD = lon) %>% 
  # add on the simulated detections
  full_join(random_detections) %>% 
  filter(!is.na(detection)) %>% 
  # convert everything to radians
  mutate(detection = activity::solartime(detection, YCOORD, XCOORD, tz = 6, 
                              format = "%Y-%m-%d %H:%M:%S")$solar,
         nightEnd = activity::solartime(nightEnd, YCOORD, XCOORD, tz = 6, 
                              format = "%Y-%m-%d %H:%M:%S")$solar,
         goldenHourEnd = activity::solartime(goldenHourEnd, YCOORD, XCOORD, tz = 6, 
                                   format = "%Y-%m-%d %H:%M:%S")$solar,
         goldenHour = activity::solartime(goldenHour, YCOORD, XCOORD, tz = 6, 
                                format = "%Y-%m-%d %H:%M:%S")$solar,
         night = activity::solartime(night, YCOORD, XCOORD, tz = 6, 
                           format = "%Y-%m-%d %H:%M:%S")$solar)

# empty list for results
period_proportions <- list()
for(i in 1:(length(dates) - 1)){
  
  current_date <- dates[i]
  
  # fit activity model to simulated detection data
  act <- activity::fitact(filter(suntimes,
                                 date == current_date)
                          %>% pull(detection))
  
  # grab just the probability density function of the activity model
  pdf <- act@pdf
  # points at which PDF is evaluated
  xx <- pdf[,1]
  # binwdith
  dx <- xx[2L] - xx[1L]
  # probability density for each bin
  yy <- pdf[,2]
  
  # get the mean period endpoints for points, since activity model
  # pools detection data from all points
  period_times <- filter(suntimes, date == current_date) %>% 
    dplyr::select(nightEnd:night) %>% 
    summarise_all(mean) 
  
  # proportion dawn-active
  # dawn is the time between night end and end of morning "golden hour"
  dawn <- sum(yy[(xx > period_times$nightEnd & xx < period_times$goldenHourEnd)])*dx
  # proportion dusk-active
  # dusk is time between start of evening "golden hour" and night start
  dusk <- sum(yy[(xx > period_times$goldenHour & xx < period_times$night)])*dx
  #proportion diurnal
  #day is between end of morning golden hour and start of evening golden hour
  diurnal <- sum(yy[(xx > period_times$goldenHourEnd & xx < period_times$goldenHour)])*dx
  #proportion nocturnal
  # night is between midnight and night end plus nigth start until midnight
  nocturnal <- sum(yy[(xx > 0 & xx < period_times$nightEnd) | (xx > period_times$night)])*dx

  period_proportions[i] <- list(tibble(
    date = current_date, 
    pdawn = dawn, 
    pdusk = dusk, 
    pdiurnal = diurnal, 
    pnocturnal = nocturnal))
}

# unpack list into a dataframe
# each day of survey is a row
# columns for proportions of activity within dawn, dusk, day, and night periods
( period_proportions_df <- bind_rows(period_proportions) )
