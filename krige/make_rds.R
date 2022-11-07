
# Imports
library('geoR')

model <- 'model4.6.7'

k0_dat_csv <- read.csv(file = sprintf("/Users/tnye/kappa/krige/%s_culled_coords.csv", model))

lon <- k0_dat_csv['Longitude']
lat <- k0_dat_csv['Latitude']
UTMx <- k0_dat_csv['UTMx']
UTMy <- k0_dat_csv['UTMy']
log_k <- k0_dat_csv['log10_kappa']
kappa <- k0_dat_csv['kappa']
covariates <- numeric(length(log_k))


### Longitude Latitude
df <- data.frame("Longitude"=lon, "Latitude"=lat, "data"=kappa, "covariates"=covariates)

data <- as.geodata(df, coords.col=1:2, data.col=3, covar.col=4,
                   "station_names"=df$station_names, data.names=NULL, covar.names="BAY")

saveRDS(data, file = sprintf("/Users/tnye/kappa/krige/%s_culled_my_data_lonlat.rds", model))
  

### UTM
df <- data.frame("UTMx"=UTMx, "UTMy"=UTMy, "data"=kappa, "covariates"=covariates)

data <- as.geodata(df, coords.col=1:2, data.col=3, covar.col=4,
                   "station_names"=df$station_names, data.names=NULL, covar.names="BAY")

saveRDS(data, file = sprintf("/Users/tnye/kappa/krige/%s_culled_my_data_utm.rds", model))
