rm(list = ls())

rb <- read.csv("LiberiaRemoData.csv")

# log(p(x)/(1+p(x))) = intercept + reg.coef * log(elevation) + S(x) + Z
# S(x) to be a sptial Gaussian process with exponential correlation function
# cov(S(x), S(x')) = sigma^2*exp(-||x-x'||/phi)

# Prameters to guess
# intercept and reg. coefficient
# variance (sigma^2) and scale the spatial correlation (phi)
# variance (tau^2) of the nugget effect Z

beta.guess <- coef(glm(cbind(npos,ntest-npos) ~ log(elevation),
                  family=binomial,data=rb))

spat.corr.diagnostic(npos~log(elevation),
                     units.m=~ntest,data=rb,
                     coords=~I(utm_x/1000)+I(utm_y/1000),
                     likelihood = "Binomial",
                     lse.variogram = TRUE)

sigma2.guess <- 0.02
phi.guess <- 36.95
tau2.guess <-0.01

#?binomial.logistic.MCML

# Vector of guessed parameters
# NOTE: The order must be respected
# The order is this: regression coefficients, variance of the Gaussian process, 
# scale of the spatial correlation, 
# variance of the nugget
par0.rb <- c(beta.guess , sigma2.guess , phi.guess , tau2.guess)

#?control.mcmc.MCML
mcml <- control.mcmc.MCML(n.sim=10000,burnin=2000,thin=8)

# par0 corresponds to theta0 on the slides
fit.mle <-
binomial.logistic.MCML(npos ~ log(elevation),
                       units.m = ~ ntest,
                       coords=~I(utm_x/1000)+I(utm_y/1000),
                       par0=par0.rb,control.mcmc = mcml,
                       kappa=0.5,
                       start.cov.pars = c(phi.guess,tau2.guess/sigma2.guess),
                       data=rb,method="nlminb")

summary(fit.mle,log.cov.pars=FALSE)

# log(p(x)/(1+p(x))) = intercept + reg.coef * log(elevation) + S(x) 
par0.rb <- c(beta.guess,sigma2.guess,phi.guess)
fit.mle.no.nugget <-
  binomial.logistic.MCML(npos ~ I(log(elevation)),
                         units.m = ~ ntest,
                         coords=~I(utm_x/1000)+I(utm_y/1000),
                         par0=par0.rb,control.mcmc = mcml,
                         kappa=0.5,
                         fixed.rel.nugget = 0,
                         start.cov.pars = phi.guess,
                         data=rb,method="nlminb")

par0.rb <- coef(fit.mle.no.nugget)
phi.guess <- par0.rb["phi"]
fit.mle.no.nugget <-
  binomial.logistic.MCML(npos ~ I(log(elevation)),
                         units.m = ~ ntest,
                         coords=~I(utm_x/1000)+I(utm_y/1000),
                         par0=par0.rb,control.mcmc = mcml,
                         kappa=0.5,
                         fixed.rel.nugget = 0,
                         start.cov.pars = phi.guess,
                         data=rb,method="nlminb")


liberia.bndrs <- read.csv("Liberia_bndrs.csv")
library(splancs)  
liberia.grid <- gridpts(
                as.matrix(liberia.bndrs[,c("utm_x","utm_y")])/1000,
                        xs=3,ys=3)  

elevation <- raster("./LBR_alt/LBR_alt.gri")
elevation <- projectRaster(elevation,crs=CRS("+init=epsg:32629"))
elevation.pred <- extract(elevation,liberia.grid*1000)

ind.na <- which(is.na(elevation.pred))
liberia.grid <- liberia.grid[-ind.na,]
elevation.pred <- elevation.pred[-ind.na]

predictors.rb <- data.frame(elevation=elevation.pred)
  
pred.mle <- 
spatial.pred.binomial.MCML(fit.mle,grid.pred = liberia.grid,
                           predictors = predictors.rb,
                           control.mcmc = mcml,
                           scale.predictions = c(
                             "logit","prevalence"),
                           thresholds = 0.2,
                           scale.thresholds = "prevalence")

plot(pred.mle,"prevalence","predictions")
plot(pred.mle,summary="exceedance.prob")

plot(pred.mle,"prevalence","predictions", col = colorRampPalette(c("dark green", "green", "yellow", "orange", "red"))(10))
