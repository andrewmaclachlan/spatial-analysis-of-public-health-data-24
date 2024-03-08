rm(list=ls())

library(PrevMap)
data("loaloa")

loaloa$e.logit <- log((loaloa$NO_INF+0.5)/(loaloa$NO_EXAM-loaloa$NO_INF+0.5))

library(ggplot2)

elev.set <- c(600,650,750,800)

glm600 <- glm(cbind(NO_INF,NO_EXAM-NO_INF) ~ ELEVATION + 
                I((ELEVATION-650)*(ELEVATION>650)),
              data=loaloa,family = binomial)

glm700 <- glm(cbind(NO_INF,NO_EXAM-NO_INF) ~ ELEVATION + 
                I((ELEVATION-700)*(ELEVATION>700)),
              data=loaloa,family = binomial)

glm800 <- glm(cbind(NO_INF,NO_EXAM-NO_INF) ~ ELEVATION + 
                I((ELEVATION-750)*(ELEVATION>750)),
              data=loaloa,family = binomial)

logLik(glm600)
logLik(glm700)
logLik(glm800)




glm.quad <- glm(cbind(NO_INF,NO_EXAM-NO_INF) ~ ELEVATION + 
             I(ELEVATION^2),
           data=loaloa,family = binomial)
summary(glm.quad)
beta.hat <- coef(glm.quad)
-beta.hat[2]/(2*beta.hat[3])



plot.elev <- ggplot(loaloa, aes(x = ELEVATION, 
                             y = e.logit)) + geom_point() +
  labs(x="Elevation",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x + I(x^2),
              col="red",lty="dashed",se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x + I((x-700)*(x>700)),
              col="green",lty="dashed",se=FALSE)


plot.ndvi <- ggplot(loaloa, aes(x = MEAN9901, 
                                y = e.logit)) + geom_point() +
  labs(x="Mean NDVI (1999-2001)",y="Empirical logit")+
  stat_smooth(method = "gam", formula = y ~ s(x),se=FALSE)+
  stat_smooth(method = "lm", formula = y ~ x ,
              col="red",lty="dashed",se=FALSE)
  

library(lme4)
loaloa$ID <- 1:nrow(loaloa)
glmer.fit <- glmer(cbind(NO_INF,NO_EXAM-NO_INF) ~ 
                   MEAN9901+(1|ID),
                   data = loaloa, family = binomial,
                   nAGQ = 100)
summary(glmer.fit)

spat.corr.diagnostic(NO_INF ~ MEAN9901, 
                     units.m = ~ NO_EXAM,
                     coords = ~ LONGITUDE + LATITUDE,
                     data = loaloa,
                     uvec = seq(0.01, 3.5, length = 15),
                     likelihood = "Binomial",
                     which.test = "variogram")


