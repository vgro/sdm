# -------------------------------------------------
# script to run SDM without land use change
#------------------------------------------------

wd <- setwd("Z:/Modelcomparison/00_data")
wd

library(dismo)
library(rgdal)  
library(raster) 
library(sp)
library(maptools)

#------------------------------------------------
# import data
#------------------------------------------------

# read in predictors
files <- list.files(pattern='tif')
files

predictors <- stack(files)
predictors
names(predictors)
plot(predictors,xlim=c(-20,50),ylim=c(-40,20))
pairs(predictors)

# read in species occurrence
data(wrld_simpl)
species_occurrence <- read.csv("Z:/Modelcomparison/occurrence/bushmeat/occurence_ready/Concern/Acinonyx_jubatus.csv", sep=',')
coordinates(species_occurrence) <- ~lon+lat
crs(species_occurrence) <- "+proj=longlat +datum=WGS84 +no_defs"
species_occurrence

# plot occurrence on map of first predictor
plot(predictors, 1)
plot(wrld_simpl, add=TRUE)
points(species_occurrence, col='blue',pch=19)

# set species extent
ext <- extent(species_occurrence)
ext

#------------------------------------------------------------------
# create data subsets
#------------------------------------------------------------------

# occurence training and test set
occ=cbind.data.frame(species_occurrence$lon,species_occurrence$lat) #first, just make a data frame of latitudes and longitudes for the model
group <- kfold(occ, k=5) # add an index that makes five random groups of observations
pres_test <- occ[group == 1, ] # hold out one fifth as test data
pres_train <- occ[group != 1, ] # the other four fifths are training data
colnames(pres_train) = c('lon', 'lat')
colnames(pres_test) = c('lon', 'lat')

# Background training and test set
set.seed(0)
backgr <- randomPoints(predictors, 500, extf=1.25) # set background points
group <- kfold(backgr, 5)
backg_train <- backgr[group != 1, ]
backg_test <- backgr[group == 1, ]
colnames(backg_train) = c('lon', 'lat')
colnames(backg_test) = c('lon', 'lat')

# plot presence and absence values
r <- raster(species_occurrence, 1)
data(wrld_simpl)
plot(wrld_simpl, add=TRUE)
plot(!is.na(r), col=c('white', 'light grey'), legend=FALSE)
plot(ext, add=TRUE, col='red', lwd=2)
points(backg_train, pch='-', cex=1, col='yellow')
points(backg_test, pch='-',  cex=1, col='black')
points(pres_train, pch= '+', col='green')
points(pres_test, pch='+', col='blue')

#------------------------------------------------------------------
# model ensemble:
# BIOCLIM
# Generalized linear regression
# MAXENT
# random forest
# boosted regression trees
#------------------------------------------------------------------

# BIOCLIM
# run model
bc <- bioclim(predictors, pres_train) #plot(bc, a=1, b=2, p=0.85)

# evaluate model
e <- evaluate(pres_test, backg_test, bc, predictors)
e

#Find a threshold, 'spec_sens'= the sum of the sensitivity (true positive rate) and specificity (true negative rate) is highest 
# see also https://rdrr.io/cran/dismo/man/threshold.html for alternatives
tr_bc <- threshold(e, 'spec_sens')  
tr_bc

# predict suitability
predict_bioclim <- predict(predictors, bc, ext=ext, progress='')
predict_bioclim

# plot raw values and presence/absence
par(mfrow=c(1,2))
plot(predict_bioclim, main='Bioclim, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(predict_bioclim > tr_bc, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(species_occurrence$lon, species_occurrence$lat, pch='+')


# GENERALIZED LINEAR MODEL

# needs presence and absence, so need to extract values
train <- rbind(pres_train, backg_train)
set_frame_train <- c(rep(1, nrow(pres_train)), rep(0, nrow(backg_train)))
envtrain <- extract(predictors, train)
envtrain <- data.frame( cbind(pa=set_frame_train, envtrain) )
head(envtrain)
testpres <- data.frame( extract(predictors, pres_test) )
testbackg <- data.frame( extract(predictors, backg_test) )

# logistic regression:
gm1 <- glm(pa ~ CHELSA_annualprecip1 + CHELSA_isothermaly1 + CHELSA_meandiurnalrange1 + CHELSA_meantemp1 + CHELSA_precipseasonality1
           + CHELSA_tempannualrange1 + CHELSA_tempmaxmax1 + CHELSA_tempminmin1 + CHELSA_tempseasonality1 + EarhEnv_elevation
           + EarthEnv_eastness + EarthEnv_northness + EarthEnv_slope + GSW_wateroccurence,
           family = binomial(link = "logit"), data=envtrain)

summary(gm1)
coef(gm1)

# correlation between predictors
library(tidyverse)
library(car)
car::vif(gm1)

# simplified model MANUAL SO FAR
gm2 <- glm(pa ~ CHELSA_annualprecip1 + CHELSA_precipseasonality1+ EarhEnv_elevation + EarthEnv_eastness 
           + EarthEnv_northness + EarthEnv_slope + GSW_wateroccurence,
           family = binomial(link = "logit"), data=envtrain)
summary(gm2)
coef(gm2)
tr_glm <- threshold(ge2, 'spec_sens')

# evaluate models
ge1 <- evaluate(testpres, testbackg, gm1)
ge2 <- evaluate(testpres, testbackg, gm2)
ge1
ge2

# predict suitability
predict_glm <- predict(predictors, gm2, ext=ext)

# plot raw values and 
par(mfrow=c(1,2))
plot(predict_glm, main='GLM/gaussian, raw values')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(predict_glm > tr_glm, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(species_occurrence$lon, species_occurrence$lat, pch='+', cex=1)
points(backg_train, pch='-', cex=1)


# MAXENT
jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')

# run model
ma <- maxent(predictors, pres_train) #note we just using the training data
ma

plot(ma) # relative importance of the different predictors in the final model
response(ma) #how does the likelihood of species occurrence respond to variation in these climatic conditions

# predict favorable habitats
predict_maxent <- predict(ma, predictors,ext=ext)

# evaluate model
e1 <- evaluate(ma, p=pres_test, a=backgr, x=predictors) #generate diagnostics and the AUC, using our test data as our presences against our pseudoabsences.
plot(e1, 'ROC')

# plot raw values and presence
plot(predict_maxent, main="Predicted Suitability")
map('worldHires', fill=FALSE, add=TRUE)
points(species_occurrence$lon, species_occurrence$lat, pch=19, cex=0.2)


# RANDOM FOREST
library(randomForest)

model <- (pa ~ CHELSA_annualprecip1 + CHELSA_isothermaly1 + CHELSA_meandiurnalrange1 + CHELSA_meantemp1 + CHELSA_precipseasonality1
          + CHELSA_tempannualrange1 + CHELSA_tempmaxmax1 + CHELSA_tempminmin1 + CHELSA_tempseasonality1 + EarhEnv_elevation
          + EarthEnv_eastness + EarthEnv_northness + EarthEnv_slope + GSW_wateroccurence)

rf2 <- randomForest(model, data=envtrain,na.action=na.roughfix)

erf <- evaluate(testpres, testbackg, rf2)
erf

predict_rf <- predict(predictors, rf2, ext=ext)
tr_rf <- threshold(erf, 'spec_sens')

par(mfrow=c(1,2))
plot(predict_rf, main='Random Forest, regression')
plot(wrld_simpl, add=TRUE, border='dark grey')
plot(predict_rf > tr_rf, main='presence/absence')
plot(wrld_simpl, add=TRUE, border='dark grey')
points(pres_train, pch='+')
points(backg_train, pch='-', cex=0.25)


# BOOSTED REGRESSION TREES
library(gbm)

brt_model <- gbm.step(data=envtrain, gbm.x = 2:10, gbm.y = 1) #Build initial model
predict_brt <- predict(predictors, brt_model, n.trees=brt_model$gbm.call$best.trees, type="response",ext=ext) #Make predictions across raster layer
plot(predict_brt, main="BRT prediction - full") #Make plot

ebrt<-  evaluate()

# MAHALANOBIS DISTANCE
#mm <- mahal(predictors, pres_train)
#e <- evaluate(pres_test, backg_test, mm, predictors)
#e

#predict_mahalanobis = predict(predictors, mm, ext=ext, progress='')
#par(mfrow=c(1,2))
#predict_mahalanobis[predict_mahalanobis < -10] <- -10
#plot(predict_mahalanobis, main='Mahalanobis distance')
#plot(wrld_simpl, add=TRUE, border='dark grey')
#tr <- threshold(e, 'spec_sens')
#plot(predict_mahalanobis > tr, main='presence/absence')
#plot(wrld_simpl, add=TRUE, border='dark grey')
#points(pres_train, pch='+')



#--------------------------------------------------------------------------------------
# COMBINING MODEL PREDICTIONS
#--------------------------------------------------------------------------------------

predict_bioclim
predict_glm
predict_maxent
predict_rf

models <- stack(predict_bioclim, predict_glm, predict_maxent, predict_rf,predict_brt)
names(models) <- c("bioclim", "glm", "maxent", "random forest", "boosted regression tree")
plot(models)


auc <- sapply(list(e,ge2,e1,erf,ebrt), function(x) x@auc)
w <- (auc-0.5)^2
m2 <- weighted.mean( models, w)
plot(m2, main='weighted mean of 5 models')
points(species_occurrence$lon, species_occurrence$lat, pch=19, cex=0.2)







