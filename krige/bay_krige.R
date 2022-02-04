##################################################################################
# R code to read station k0 data from Van Houtte et al. (2018), "A continuous 
# near-surface S-wave attenuation map of New Zealand", Geophysical Journal 
# International, https://doi.org/10.1093/gji/ggx559, and calculates the continuous 
# k0 maps using kriging. 
# Please refer to the article, for explanations of 'Model 1' and 'Model 2'
#
# Requires R packages {geoR} and {pracma}.
#
##################################################################################

# Imports
library('geoR')
library('pracma')

model <- 'model2'
nug.model <- 0.0043

# Read data:
# - Coordinates are NZTM northings and eastings, converted to km for convenience.
# - Data are log10(k0)
k0_dat <- readRDS(sprintf("/Users/tnye/kappa/krige/%s_my_data.rds", model))
grid <- read.csv("/Users/tnye/kappa/krige/grid_lonlat.txt")

# Initial values for regression
ini.nug <- 0
ini.phi <- 100
ini.sill <- 0.1

# Matern order (fixed)
theta <- 0.5 

# Model 1 - solving for nugget

model <- likfit(k0_dat, cov.model="matern", ini.cov.pars = c(ini.sill, ini.phi), 
                 nug = ini.nug, fix.nugget=F, kappa=theta, fix.kappa=T,trend="cte")

# Prediction grid in UTM coordinates, i.e. northing and easting, converted to km
# pred.grid<- expand.grid(seq(500000, 674000, l=175), 
#                         seq(4028316, 4262316, l=235)) / 1000

#pred.grid<- expand.grid(seq(-123, -121, l=201), 
#                        seq(36.4, 38.5, l=211))

# Read digitised 'whole TVZ' model of Wilson et al. (1995), coordinates in km
#whole.tvz=read.table("/Users/tnye/code/kappa/codes/whole_tvz_ne.txt", header = TRUE)
# whole.bay=read.table("/Users/tnye/kappa/krige/ne.txt", header = TRUE)
# fbay<-ifelse(inpolygon(pred.grid$Var1, pred.grid$Var2, 
#                        whole.bay$Easting/1000, whole.bay$Northing/1000), 1, 0)
whole.bay=read.table("/Users/tnye/kappa/krige/lonlat.txt", header = TRUE)
# fbay<-ifelse(inpolygon(pred.grid$Var1, pred.grid$Var2, 
#                        whole.bay$Longitude, whole.bay$Latitude), 1, 0)
fbay<-ifelse(inpolygon(grid$Longitude, grid$Latitude, 
                       whole.bay$Longitude, whole.bay$Latitude), 1, 0)
# Prediction
# nug.model <- model$tausq
#nug.model <- 0.0026
theta.model <- 0.5
sill.model <- model$sigmasq
phi.model <- model$phi
psiA.model <- model$aniso.pars[1]
psiR.model <- model$aniso.pars[2]

kc.model <- krige.conv(k0_dat, loc = grid,
                       krige = krige.control(type.krige="ok", cov.model="matern",
                                             cov.pars=c(sill.model, phi.model),
                                             kappa=theta.model, nugget=nug.model,
                                             aniso.pars=c(psiA.model, psiR.model),
                                             trend.d="cte", trend.l="cte"),
                       output = output.control(signal=T))

# Units of grid back to m
Nor <- kc.model$Var2 * 1000
Eas <- kc.model$Var1 * 1000
k0.model <- 10^kc.model$predict
logk0.stddev <- sqrt(kc.model$krige.var + nug.model)

df <- data.frame("Latitude"=grid[,2], "Longitude"=grid[,1], "pred_k0"=k0.model, "log10k0_stddev"=logk0.stddev)
write.csv(df, sprintf("/Users/tnye/kappa/krige/%s_krige_k0_lonlat.csv", model), row.names = FALSE)
#(df, "/Users/tnye/kappa/krige/model2_krige_k0_lonlat.csv", row.names = FALSE)



