##################################################################################
# R code to read station k0 data and produce a spatially-continuous map for the 
# San Francisco Bay area. Code modified from 
# https://github.com/cvanhoutte-zz/kappa/blob/master/codes/krige.R
#
# Requires R packages {geoR} and {pracma}.
#
##################################################################################

# Imports
library('geoR')
library('pracma')
library('gstat')
library(sp)

# Name of model
model_name <- 'model4.6.7'

# Read data for kappa estimates from rds file and turn into a dataframe:
# - Coordinates are degrees lon and lat
k0rds <- readRDS(sprintf("/Users/tnye/kappa/krige/model4.6.7/%s_culled_my_data_lonlat.rds", model_name))
kappa <- k0rds$data
longitude <- k0rds$coords[,1]
latitude <- k0rds$coords[,2]
k0.data <- data.frame(longitude, latitude, kappa)

# Define coordinates
coordinates(k0.data) = ~longitude+latitude

# Initial values for regression (from Van Houtte et al., 2018)
ini.nug <- 0
ini.sill <- 0.0002
ini.range <- 0.4

# Matern order (fixed) (from Van Houtte et al., 2018)
theta <- 0.5 

# Make and plot variogram
k0.vgm = variogram(kappa~1, k0.data)
k0.fit = fit.variogram(k0.vgm, model = vgm(ini.sill, 'Sph', ini.range, nugget=ini.nug)) 

# Save plot of variogram
png(file=sprintf('/Users/tnye/kappa/plots/paper/k0_semivariogram_%s_culled.png',model_name),
    width=600, height=350, res=300)
plot(k0.vgm, k0.fit, xlab="Distance (deg)", ylab="Semivariance",cex.lab=.5,cex.axis=.5)
dev.off()

png(file=sprintf('/Users/tnye/kappa/plots/paper/k0_semivariogram_%s_culled.png',model_name),
    width=600, height=350, res=300)
plot(k0.vgm, k0.fit, xlab="Distance (deg)", ylab="Semivariance",cex.lab=.1,cex.axis=.1)
dev.off()

# Kriging 1 --------------------------------------------------------------------
grid <- read.csv("/Users/tnye/kappa/krige/grid_lonlat2.txt")
coordinates(grid) = ~Longitude+Latitude
k0.kriged = krige(kappa~1, k0.data, grid, model = k0.fit)
spplot(k0.kriged["var1.pred"])

# Krige Predictions (in linear, not log10 units)
k0.pred <- k0.kriged$var1.pred

# Krige Standard deviations
k0.stddev <- sqrt(k0.kriged$var1.var)

# Make Dataframe
df <- data.frame("Latitude"=grid$Latitude, "Longitude"=grid$Longitude, "pred_k0"=k0.pred, "k0_stddev"=k0.stddev)
write.csv(df, sprintf("/Users/tnye/kappa/krige/%s/%s_culled_krige_k0_linear.csv", model_name, model_name), row.names = FALSE)


