
# Imports
library('geoR')

k0_dat_csv <- read.csv(file = '/Users/tnye/kappa/krige/model1_coords.csv')

# coords <- k0_dat_csv[,3:4]
# log_k <- k0_dat_csv[,5]
# covariates <- numeric(length(log_k))
# kappa <- 0.5
# model <- "matern"
# stations <- k0_dat_csv[,2]

# coords <- k0_dat_csv[,2:3]
coords <- k0_dat_csv[,4:5]
log_k <- k0_dat_csv[,6]
covariates <- numeric(length(log_k))

# df <- data.frame("UTMx"=coords[,1], "UTMy"=coords[,2], "data"=log_k, "covariates"=covariates)
df <- data.frame("Longitude"=coords[,1], "Latitude"=coords[,2], "data"=log_k, "covariates"=covariates)

# data <- as.geodata(df, coords.col=1:2, data.col=3, covariates.col=4,
#                    "station_names"=df$station_names, data.names=NULL, covar.col=NULL,
#                    covar.names="obj.names")

data <- as.geodata(df, coords.col=1:2, data.col=3, covar.col=4,
                   "station_names"=df$station_names, data.names=NULL, covar.names="BAY")

saveRDS(data, file = "/Users/tnye/kappa/krige/model1_my_data.rds")
