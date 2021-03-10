
#  Malawi, grain nutrient concentrations modelled with (i) environmental covariates
#  or (ii) soil properties.
#
#  The analyses on the Malawi data are done using the geoR package.  This is 
#  because the data set is much larger than in Ethiopia, which is computationally
#  demanding, and also, because Malawi is narrow from East to West the coordinates
#  can be projected onto the Universal Transverse Mercator projection for a single
#  zone.  This allows the computations to be done more rapidly.
#
#  This script is for maize Se concentration, as with the Ethiopian analysis
#  each variable may have a different sequence of models fitted, depending
#  on whether the predictors are accepted or rejected at each stage.
#
#############################################################################

library(xlsx)
library(maps)
library(mapdata) 
library(sp)
library(geoR)
source("CEPHaStat_2.R")

#############################################################################

data.df<-read.xlsx("Malawi_grain_soil.xlsx",sheetIndex=1,header=T)
data.df<-data.df[which(data.df$Crop=="Maize"),]
data.df<-data.df[complete.cases(data.df$Se_triplequad),]
data.df<-data.df[complete.cases(data.df$pH_Wa),]
data.df<-data.df[complete.cases(data.df$Org_pct),]


##########################################################################
#  
#  Analyses with soil properties
#

#  Exploratory analysis

par(mfrow=c(1,2))
par(mar=c(5,5,3,3))
plot(data.df$pH_Wa,data.df$Se_triplequad,pch=16,xlab="Soil pH",
ylab=expression("Maize grain Se/ mg kg"^{-1}),log="y")

par(mar=c(5,5,3,3))
plot(data.df$Org_pct,data.df$Se_triplequad,pch=16,xlab="Soil organic carbon/ %",
ylab=expression("Maize grain Se/ mg kg"^{-1}),log="xy")

mod<-lm(log(Se_triplequad)~pH_Wa+log(Org_pct),data=data.df)

summa(mod$residuals)
summaplot(mod$residuals)

par(mfrow=c(2,2))
par(mar=c(5,5,2,2))
plot(data.df$BIO12,data.df$Se_triplequad,pch=16,xlab="Mean annual precipitation /mm",
ylab=expression("Maize grain Se/ mg kg"^{-1}),log="y")

par(mar=c(5,5,2,2))
plot((data.df$BIO1)/10,data.df$Se_triplequad,pch=16,xlab="Mean annual temperature /°C",
ylab=expression("Maize grain Se/ mg kg"^{-1}),log="y")

par(mar=c(5,5,2,2))
plot((data.df$TIM),data.df$Se_triplequad,pch=16,xlab="Topographic Index",
ylab=expression("Maize grain Se/ mg kg"^{-1}),log="y")


mod<-lm(log(Se_triplequad)~(TIM),data=data.df)

summary(mod)
anova(mod)
summa(mod$residuals)
summaplot(mod$residuals)


###########################################################################
#
#  Project the data onto the UTM zone 36S
#
UTM36S="+init=epsg:32736"

locs<-cbind(data.df$Longitude,data.df$Latitude)
loc_sp<-SpatialPoints(locs, proj4string=CRS("+proj=longlat +datum=WGS84"))
loc_UTM<- spTransform(loc_sp, CRS(UTM36S))


#  Create a dataframe with the UTM coordinates, grain nutrient conc.
#  and soil properties

gdata<-data.frame(cbind(loc_UTM@coords,data.df$Se_triplequad,
data.df$pH_Wa,data.df$Org_pct))

#
#  Next, make a geodata object for the geoR library procedures
#

geod<-as.geodata(gdata,covar.col = c(4,5))

#
#  Fit a null model in which the only fixed effect is a constant mean
#

mod0<-likfit(geod,
trend="cte",cov.model="exponential",ini.cov.pars=c(1,30.0),lambda=0,
lik.method="ML")

#Extract the maximized log-likelihood for the null model

l0<-mod0$loglik  #2701.842

# Next, include the first soil property (pH) as a proposed predictor
mod1<-likfit(geod,
trend=~V4,cov.model="exponential",ini.cov.pars=c(1,30.0),lambda=0,
lik.method="ML")

#Extract the maximized log-likelihood

l1<-mod1$loglik  #2709.729

# Compute the log-likelihood ratio statistic for adding pH, and the p-value
# for the null hypothesis,

L1<-2*(l1-l0)  #15.77  
pval<-1-pchisq(L1,1) #p=7.1 x 10^-5

#Extract the fixed effects coefficients
mod1$beta
#  intercept     covar1 
#-4.7207352  0.1895304 

#Extract the standard error of the coefficient for pH

sqrt(mod1$beta.var[2,2]) #0.0474

#  Next add (log) of soil organic carbon as a potential covariate

mod2<-likfit(geod,
trend=~V4+log(V5),cov.model="exponential",ini.cov.pars=c(1,30.0),lambda=0,
lik.method="ML")

#Extract the maximized log-likelihood

l2<-mod2$loglik  #2710.85

# Compute the log-likelihood ratio statistic for adding organic carbon, and the p-value
# for the null hypothesis,

L2<-2*(l2-l1)  #2.24  
pval<-1-pchisq(L2,1) #p=0.134

###################################################################################
#
#  Environmental covariates

#Create a dataframe with the grain nutrient concentration and environmental covariates

gdata<-data.frame(cbind(loc_UTM@coords,data.df$Se_triplequad,
data.df$BIO12,data.df$BIO1,data.df$TIM))

geod<-as.geodata(gdata,covar.col = c(4,5,6))

# Fit the null model 

emod0<-likfit(geod,trend="cte",cov.model="exponential",ini.cov.pars=c(1,30.0),lambda=0,
lik.method="ML")

# Extract the maximized log-likelihood for the null model

el0<-emod0$loglik  #2693.729


# Add the first environmental covariate (Downscaled mean annual precipitation, MAP)

emod1<-likfit(geod,
trend=~V4,cov.model="exponential",ini.cov.pars=c(1,30.0),lambda=0,
lik.method="ML")

# Extract the maximized log-likelihood 

el1<-emod1$loglik  #2693.74

# Compute the log-likelihood ratio statistic for adding MAP, and the p-value
# for the null hypothesis,

Le1<-2*(el1-el0)  #0.02  
pval<-1-pchisq(Le1,1) #p=0.88

#  There is no evidence to retain MAP in the model, so the next model as Downscaled mean 
#  annual temperature (MAT) as the only covariate

emod2<-likfit(geod,
trend=~V5,cov.model="exponential",ini.cov.pars=c(1,30.0),lambda=0,
lik.method="ML")

# Extract the maximized log-likelihood 

el2<-emod2$loglik  #2701.375

# Compute the log-likelihood ratio statistic for adding MAP, and the p-value
# for the null hypothesis,

Le2<-2*(el2-el0)  #15.29  
pval<-1-pchisq(Le2,1) #9.209108e-05

# We have evidence to retain MAP as a predictor

#Extract the model parameters and the standard error of the effect for MAT

emod2$beta
#  intercept      covar1 
#-6.36022190  0.01321597 

#se of parameter
sqrt(emod2$beta.var[2,2]) #0.00320

# Add Terrain Index (TIM) as a predictor

emod3<-likfit(geod,
trend=~V5+V6,cov.model="exponential",ini.cov.pars=c(1,30.0),lambda=0,
lik.method="ML")

# Extract the maximized log-likelihood 

el3<-emod3$loglik  # 2701.681

# Compute the log-likelihood ratio statistic for adding MAP, and the p-value
# for the null hypothesis,

Le3<-2*(el3-el2)  #0.61  
pval<-1-pchisq(Le3,1) #0.433

#  There is no evidence to retain TIM in the model.
