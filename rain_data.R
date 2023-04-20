
# get rain data ####

require(rnoaa)
require(lubridate)

lat_lon_df <- data.frame(id = "TZ000063789", 
                         lat = -3.333, 
                         lon = 36.633)

mon_near_arusha <- 
  meteo_nearby_stations(
    lat_lon_df = lat_lon_df,
    lat_colname = "lat",
    lon_colname = "lon",
    var = "PRCP",
    year_min = 1980,
    year_max = 1989,
    limit = 10,
  )

arusha_prcp_dat <- 
  meteo_pull_monitors(
    monitors = mon_near_arusha$TZ000063789$id[1],
    date_min = "1980-02-01",
    date_max = "1989-12-31",
    var = "PRCP"
  )

arusha_prcp_dat %>% 
  mutate(year = year(date),
         month = month(date)) %>%
  group_by(year, month) %>% 
  summarise(prcp = sum(prcp, na.rm = TRUE) / 10) %>% 
  ggplot(aes(factor(month), prcp))+
  geom_col()+
  facet_wrap(~year, nrow=2)+
  labs(y = "Precipitation [mm]",
       x = "Month")+
  theme_bw()+
  theme(axis.text.x = element_text(size = 10), 
        axis.title = element_text(size = 16), 
        legend.text = element_text(size = 16), 
        legend.title = element_text(size = 16), 
        strip.text = element_text(size = 12))
  
ggsave("precipitation 1980 to 89.pdf")

# from the precipitation data when will eggs hatch
precipdata <- as.data.frame(arusha_prcp_dat$prcp)
precipdata$time <- 1:nrow(precipdata)
precipdata$hatch <- NA
precipdata <- do.call("rbind", replicate(10, precipdata, simplify = FALSE))
end <- nrow(precipdata)-15
# require seven days of rain, from the 7th day for 9 days eggs can hatch as per mindy's work
for(i in 1:end){
  if(precipdata[i, "arusha_prcp_dat$prcp"] > 0 & precipdata[i+1, "arusha_prcp_dat$prcp"] >0 & 
     precipdata[i+2, "arusha_prcp_dat$prcp"] >0 & precipdata[i+3, "arusha_prcp_dat$prcp"] >0 & 
     precipdata[i+4, "arusha_prcp_dat$prcp"] >0 & precipdata[i+5, "arusha_prcp_dat$prcp"] >0 & 
     precipdata[i+6, "arusha_prcp_dat$prcp"] >0){ 
    precipdata[(i+6):(i+15), "hatch"] <- 1
  }
}

precipdata$hatch[which(is.na(precipdata$hatch))] <- 0

hatchtimes <- which(precipdata$hatch>0)

hatchtimes = hatchtimes[hatchtimes %in% t]

