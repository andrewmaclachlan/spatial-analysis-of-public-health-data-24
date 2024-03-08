rm(list = ls())
library(PrevMap)
rb <- read.csv("LiberiaRemoData.csv")

# Scaling to km (optional)
rb[,c("utm_x","utm_y")] <- rb[,c("utm_x","utm_y")]/1000
rb$logit <- log((rb$npos+0.5)/(rb$ntest-rb$npos+0.5))

### POINT 1
### Standard linear model
lm.fit <- lm(logit ~ log(elevation),data=rb)

library(sf)
liberia.adm0 <- st_read("../Data/Liberia spatial data/LBR_adm/LBR_adm0.shp")
liberia.adm0 <- st_transform(liberia.adm0,crs=32629)

liberia.grid <- st_make_grid(liberia.adm0,
                             cellsize = 2000,
                             what="centers")

liberia.inout <- st_intersects(liberia.grid,
                               liberia.adm0,
                               sparse = FALSE)
liberia.grid <- liberia.grid[liberia.inout]
liberia.grid <- st_coordinates(liberia.grid)/1000

elevation <- raster("Liberia spatial data/LBR_alt/LBR_alt.gri")
elevation <- projectRaster(elevation,crs=CRS("+init=epsg:32629"))
elevation.pred <- extract(elevation,liberia.grid*1000)

ind.na <- which(is.na(elevation.pred))
liberia.grid <- liberia.grid[-ind.na,]
elevation.pred <- elevation.pred[-ind.na]
predictors.rb <- data.frame(elevation=elevation.pred)

# From the empirical logit model
# predictions of prevalence can be obtained by using the anti-logit 
# (i.e. the inverse function of the logit)
# f(x) = exp(x)/(1+exp(x)) = 1/(1+exp(-x))
pred.lm <- 1/(1+exp(-predict(lm.fit,newdata=predictors.rb)))

r.pred.lm <- rasterFromXYZ(cbind(liberia.grid,pred.lm))

plot(r.pred.lm)


#### POINT 2
### Linear geostatistical model
spat.corr.diagnostic(logit~1,
                     data=rb,
                     coords=~utm_x+utm_y,
                     likelihood = "Gaussian",
                     lse.variogram = TRUE)

sigma2.guess <- 0.5
phi.guess <- 107
tau2.guess <-0.17


fit.mle <- linear.model.MLE(logit ~ 1,
                            coords=~utm_x+utm_y,
                            kappa=0.5,
                            start.cov.pars = c(phi.guess,tau2.guess/sigma2.guess),
                         data=rb,method="nlminb")

summary(fit.mle,log.cov.pars=FALSE)


### POINT 3
pred.mle.lm <- 
  spatial.pred.linear.MLE(fit.mle,grid.pred = liberia.grid,
                          standard.errors = TRUE,
                          scale.predictions = "prevalence",n.sim.prev = 1000,
                          thresholds = 0.2, scale.thresholds = "prevalence")

plot(pred.mle.lm,"prevalence","predictions")
plot(pred.mle.lm,"prevalence","standard.errors")
plot(pred.mle.lm,summary="exceedance.prob")

### POINT 4
#### Introducing elevation

fit.mle.elev <- linear.model.MLE(logit ~ log(elevation),
                            coords=~utm_x+utm_y,
                            kappa=0.5,
                            start.cov.pars = c(phi.guess,tau2.guess/sigma2.guess),
                            data=rb,method="nlminb")

pred.mle.lm.elev <- 
  spatial.pred.linear.MLE(fit.mle.elev,grid.pred = liberia.grid,
                          predictors = predictors.rb,
                          scale.predictions = "prevalence",
                          thresholds = 0.2,n.sim.prev = 1000,
                          scale.thresholds = "prevalence")

# thresholds=0.2 allows you to compute the predictive probability that prevalence is above 0.2 
# this is also known as the exceedance probability of a 0.2 prevalence threshold

plot(pred.mle.lm.elev,"prevalence","predictions")
plot(pred.mle.lm.elev,summary="exceedance.prob")

### POINT 5 
### Binomial geostatistical model
c.mcmc <- control.mcmc.MCML(n.sim = 10000,
                            burnin=2000,
                            thin=8)

glm.fit <- glm(cbind(npos,ntest-npos) ~ log(elevation),data=rb,family = binomial)
par0 <- c(coef(glm.fit),0.2,32)
fit.bin.elev <- binomial.logistic.MCML(npos ~ log(elevation),
                                  units.m = ~ ntest,
                                  coords=~utm_x+utm_y,
                                  kappa=0.5,control.mcmc = c.mcmc,
                                  fixed.rel.nugget = 0,
                                  par0=par0,
                                  start.cov.pars = 32,
                                  data=rb,method="nlminb")


pred.mle.bin.elev <- 
  spatial.pred.binomial.MCML(fit.bin.elev,grid.pred = liberia.grid,
                             control.mcmc = c.mcmc,
                             predictors = predictors.rb,
                          scale.predictions = "prevalence",
                          thresholds = 0.2,
                          scale.thresholds = "prevalence")
par(mfrow=c(1,2))
plot(pred.mle.lm.elev$prevalence$predictions,
     pred.mle.bin.elev$prevalence$predictions,
     xlab="Linear model (empirical logit)",
     ylab="Binomial model",pch=20)
abline(0,1,col=2,lwd=2)

plot(pred.mle.lm.elev$exceedance.prob,
     pred.mle.bin.elev$exceedance.prob,
     xlab="Linear model (empirical logit)",
     ylab="Binomial model",pch=20)
abline(0,1,col=2,lwd=2)

c.mcmc <- control.mcmc.MCML(n.sim = 10000,
                            burnin=2000,
                            thin=8)

# Including the nugget
glm.fit <- glm(cbind(npos,ntest-npos) ~ log(elevation),data=rb,family = binomial)
par0 <- c(coef(glm.fit),1,32,1)
fit.bin.elev <- binomial.logistic.MCML(npos ~ log(elevation),
                                       units.m = ~ ntest,
                                       coords=~utm_x+utm_y,
                                       kappa=0.5,control.mcmc = c.mcmc,
                                       par0=par0,
                                       start.cov.pars = c(32,1),
                                       data=rb,method="nlminb")

pred.mle.bin.elev <- 
  spatial.pred.binomial.MCML(fit.bin.elev,grid.pred = liberia.grid,
                             control.mcmc = c.mcmc,
                             predictors = predictors.rb,
                             scale.predictions = "prevalence",
                             thresholds = 0.2,
                             scale.thresholds = "prevalence")

par(mfrow=c(1,2))
plot(pred.mle.lm.elev$prevalence$predictions,
     pred.mle.bin.elev$prevalence$predictions,
     xlab="Linear model (empirical logit)",
     ylab="Binomial model",pch=20)
abline(0,1,col=2,lwd=2)

plot(pred.mle.lm.elev$exceedance.prob,
     pred.mle.bin.elev$exceedance.prob,
     xlab="Linear model (empirical logit)",
     ylab="Binomial model",pch=20)
abline(0,1,col=2,lwd=2)
