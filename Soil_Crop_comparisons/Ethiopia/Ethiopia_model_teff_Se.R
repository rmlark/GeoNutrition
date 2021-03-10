#
#  Ethiopia, grain nutrient concentrations modelled with (i) environmental covariates
#  or (ii) soil properties.
#
#  This script is for the analysis of Teff Se, selected because all environmental covariates #  are chosen in the model.  For other variables different decisions are taken on
#  inclusion of covariates, which would be made by adapting this script. For example
#  with Maize S, mean annual precipitation (MAP) is retained as a predictor, but temperature
#  is not, so the final model tested has MAP and Terrain Index as predictors (and Terrain
#  index is not retained)
#
#  Note there are function commands at the bottom of this script which must be run before
#  the analyses are attempted.
#
#############################################################################
library(xlsx)
library(geosphere)
library(Matrix)
library(maps)
library(mapdata) 
library(sp)
library(geoR)
source("CEPHaStat_2.R")
#############################################################################

datain.df<-read.xlsx("Ethiopia_data.xlsx",sheetIndex=1,header=T)

element<-"Se"

varcol<-which(colnames(datain.df)==element)

data.df<-datain.df[which(datain.df[,varcol]<77777),]

cropname<-"Teff"

data.df<-data.df[which(data.df$Crop==cropname),]
data.df<-data.df[complete.cases(data.df$Se),]
data.df<-data.df[complete.cases(data.df$Bio1),]
data.df<-data.df[complete.cases(data.df$Bio12),]
names(data.df)

##########################################################################################
#
# ENVIRONMENTAL VARIABLES
#
# Bio12:	Mean annual precipitation in mm
# Bio1:	Mean annual temperature x 10 deg C
# TIM:	Terrain index
#
#  Exploratory analysis
#
#
par(mfrow=c(2,2))
par(mar=c(5,5,3,3))
plot(data.df$Bio12,data.df$Se,pch=16,xlab="Mean annual precipitation /mm",
ylab=expression("Maize grain Se/ mg kg"^{-1}),log="y")


plot((data.df$Bio1/10),data.df$Se,pch=16,xlab="Mean annual temperature /deg C",
ylab=expression("Maize grain Se/ mg kg"^{-1}),log="y")


plot((data.df$TIM),data.df$Se,pch=16,xlab="Terrain Index",
ylab=expression("Maize grain Se/ mg kg"^{-1}),log="y")

#Exploratory analysis of the residuals from am OLS model fit

mod<-lm(log(Se)~Bio12+Bio1+TIM,data=data.df)

summary(mod)
anova(mod)
summa(mod$residuals)
summaplot(mod$residuals)


####################################################################################
####################################################################################
#
#  Modelling grain nutrient concentration as a function of environmental covariates


N<-nrow(data.df)

# Make a distance matrix, and rescale to km

D<-matrix(0,nrow=N,ncol=N)

for (i in 1:(N-1)){
print(c(i,i/(N-1)))
for (j in i:N){
D[i,j]<-distVincentySphere(c(data.df$Longitude[i],data.df$Latitude[i]),
c(data.df$Longitude[j],data.df$Latitude[j]))/1000 #
D[j,i]<-D[i,j]
}
}

z<-log(data.df$Se) #Target variable, log for Se in all crops, no transformation for Zn

kap<-0.5 #Setting the smoothness parameter of the Matern correlation function to 0.5 (exponential)

######################################

# Add in Mean Annual Precipitation, first predictor considered

# Add the covariate to the design matrix

M<-matrix(c(rep(1,N),data.df$Bio12),nrow=N,ncol=2)

# Find the maximum residual log-likelihood estimates of the parameters

oo<-optim(c(100,-1,-1),rsl)
oo2<-optim(oo$par,rsl)

# find the generalized least squares estimates of the fixed effects parameters
# and their standard error (first row is intercept).

gls(oo2$par)

# Use rslred to find the log residual likelihood for the model with the added predictor
# dropped.  This is done with the projection matrix for the full model (following 
# Welham S J and Thompson R 1997 Likelihood ratio tests for fixed
# model terms using residual maximum likelihood J. R. Stat. Soc. B 59 701–14)

redlist<-2
reduced<-optim(oo2$par,rslred)

llrt(reduced$value,oo2$value,redlist) #This prints out the log-likelihood ratio statistic, 
#			its degrees of freedom, and p-value for the null hypothesis.


# Add in Mean Annual Temperature as second covariate

M<-matrix(c(rep(1,N),data.df$Bio12,data.df$Bio1),nrow=N,ncol=3)

oo<-optim(c(100,-1,-1),rsl)
oo3<-optim(oo$par,rsl)

gls(oo2$par)

redlist<-3
reduced<-optim(oo3$par,rslred)

llrt(reduced$value,oo3$value,redlist)


# Add in Terrain index as third covariate


M<-matrix(c(rep(1,N),data.df$Bio12,data.df$Bio1,data.df$TIM),nrow=N,ncol=4)

oo<-optim(c(100,-1,-1),rsl)
oo4<-optim(oo$par,rsl)

redlist<-4
reduced<-optim(oo4$par,rslred)

llrt(reduced$value,oo4$value,redlist)

gls(oo4$par)


#########################################################################
#
# SOIL VARIABLES
#

# Read in the soil data, and extract values which correspond to crop samples

N<-nrow(data.df)

soil.df<-read.xlsx("GeoN_Ethiopia_pH_OCdata.xlsx",sheetIndex=1,header=T)

soil<-matrix(nrow=N,ncol=2)
colnames(soil)<-c("pH_Wa","Org_pct")
for (i in 1:N){
soilID<-data.df$SoilSampleID[i]
soil_indx<-which(soil.df$SoilSampleID==soilID)
if(length(soil_indx)==0){
soil[i,1]<-NA
soil[i,2]<-NA
}else{
soil[i,1]<-soil.df$pH_Wa[soil_indx]
soil[i,2]<-soil.df$Org_pct[soil_indx]
}
}
data.df$pH_Wa<-soil[,1]
data.df$Org_pct<-soil[,2]

data.df<-data.df[complete.cases(data.df$pH_Wa),]
data.df<-data.df[complete.cases(data.df$Org_pct),]
N<-nrow(data.df)


#Exploratory analysis of the residuals from am OLS model fit

par(mfrow=c(1,2))
par(mar=c(5,5,3,3))
plot(data.df$pH_Wa,data.df$Se,pch=16,xlab="Soil pH in water",
ylab=expression("Teff grain Se/ mg kg"^{-1}),log="y")


plot(data.df$Org_pct,data.df$Se,pch=16,xlab="Soil organic carbon /%",
ylab=expression("Teff grain Se/ mg kg"^{-1}),log="y")


mod<-lm(log(Se)~pH_Wa+log(Org_pct),data=data.df)

summary(mod)
anova(mod)
summa(mod$residuals)
summaplot(mod$residuals)


###############################################################################
###############################################################################
#
#  Modelling grain nutrient concentration as a function of soil properties
#

N<-nrow(data.df)

# # Make a distance matrix, and rescale to km

D<-matrix(0,nrow=N,ncol=N)

for (i in 1:(N-1)){
print(c(i,i/(N-1)))
for (j in i:N){
D[i,j]<-distVincentySphere(c(data.df$Longitude[i],data.df$Latitude[i]),
c(data.df$Longitude[j],data.df$Latitude[j]))/1000 #
D[j,i]<-D[i,j]
}
}

z<-log(data.df$Se) #Target variable, log for Se in all crops, no transformation for Zn

kap<-0.5 #Setting the smoothness parameter of the Matern correlation function to 0.5 (exponential)


#  First model, pH in water as the only predictor
#   Add the soil property to the design matrix

M<-matrix(c(rep(1,N),data.df$pH_Wa),nrow=N,ncol=2)

# Find the maximum residual log-likelihood estimates of the parameters

oo<-optim(c(100,-1,-1),rsl)
oo2<-optim(oo$par,rsl)

# find the generalized least squares estimates of the fixed effects parameters
# and their standard error (first row is intercept).

gls(oo2$par)

# Use rslred to find the log residual likelihood for the model with the added predictor
# dropped.  

redlist<-2
reduced<-optim(oo2$par,rslred)

llrt(reduced$value,oo2$value,redlist) # There is evidence to retain pH as a predictor

# Now add in soil organic carbon as a second predictor

M<-matrix(c(rep(1,N),data.df$pH_Wa,data.df$Org_pct),nrow=N,ncol=3)

oo<-optim(c(100,-1,-1),rsl)
oo2<-optim(oo$par,rsl)
oo2$value

redlist<-3
reduced<-optim(oo2$par,rslred)

llrt(reduced$value,oo2$value,redlist) # No evidence to retain organic carbon in the model



##################################################################################
###################################################################################################
#
# FUNCTIONS
#
###################################################################################################

ocskew<-function(x){
octiles<-quantile(x,c(1/8,1/2,7/8))
ocs<-((octiles[3]-octiles[2])-(octiles[2]-octiles[1]))/(octiles[3]-octiles[1])
return(ocs)
}



########################################################################
#
# Residual log-likelihood following Marchant et al. 2009 notation
#
#
rsl<-function(theta){
#print(theta)
ph<-theta[1]
c0<-exp(theta[2])
c1<-exp(theta[3])

#print(c(ph,c0,c1))

A<-matern(D,ph,kap)

V<-(c0*diag(N))+(c1*A)

Bhat<-(solve(t(M)%*%solve(V)%*%M))%*%(t(M)%*%solve(V)%*%z)

rankM<-rankMatrix(M)[1]

p1<-0.5*(N-rankM)*log(2*pi)-det(M%*%t(M))
L=chol(V)
p2<-sum(log(diag(L)))
p3<-0.5*log(det(t(M)%*%solve(V)%*%M))
res<-(z-M%*%Bhat)
p4<-0.5*t(res)%*%solve(V)%*%res

return(p2+p3+p4-p1) # negative rll


}

########################################################################
#
# Residual log-likelihood following Marchant et al. 2009 notation for
# reduced model (here, model with predictor in last column of M removed)
# but projection matrix is for the full model.  M0 is the design matrix
# for the reduced model.
#
#
rslred<-function(theta){

ph<-theta[1]
c0<-exp(theta[2])
c1<-exp(theta[3])


# redlist contains list of columns to be dropped from M
M0<-M[,-redlist]

A<-matern(D,ph,kap)

V<-(c0*diag(N))+(c1*A)

Bhat<-(solve(t(M0)%*%solve(V)%*%M0))%*%(t(M0)%*%solve(V)%*%z)

rankM<-rankMatrix(M)[1]

p1<-0.5*(N-rankM)*log(2*pi)-det(M%*%t(M))
p2<-0.5*log(det(V))
p3<-0.5*log(det(t(M)%*%solve(V)%*%M))
res<-(z-M0%*%Bhat)
p4<-0.5*t(res)%*%solve(V)%*%res

return(p2+p3+p4-p1) # negative rll


}


###########################################################################
#
# extract generalized least squares estimates and their standard error
#

gls<-function(theta){

ph<-theta[1]
c0<-exp(theta[2])
c1<-exp(theta[3])

A<-matern(D,ph,kap)

V<-(c0*diag(N))+(c1*A)

Bhat<-(solve(t(M)%*%solve(V)%*%M))%*%(t(M)%*%solve(V)%*%z)

C<-solve(t(M)%*%solve(V)%*%M)

op<-(cbind(Bhat,sqrt(diag(C))))
colnames(op)<-c("Estimate","Standard error")

return(op)
}

###########################################################################
#
# extract residuals
#

resids<-function(theta){

ph<-theta[1]
c0<-theta[2]
c1<-theta[3]

A<-matern(D,ph,kap)

V<-(c0*diag(N))+(c1*A)

Bhat<-(solve(t(M)%*%solve(V)%*%M))%*%(t(M)%*%solve(V)%*%z)

res<-(z-M%*%Bhat)

return(res)
}

###########################################################################
#
# extract fitted values
#

fitval<-function(theta){

ph<-theta[1]
c0<-theta[2]
c1<-theta[3]

A<-matern(D,ph,kap)

V<-(c0*diag(N))+(c1*A)

Bhat<-(solve(t(M)%*%solve(V)%*%M))%*%(t(M)%*%solve(V)%*%z)

fitv<-(M%*%Bhat)

return(fitv)
}

###########################################################################
#
# LLRT
#

llrt<-function(nll1,nll2,redlist){

L<-2*(nll1-nll2)
df<-length(redlist)
pval<-1-pchisq(L,df)

return(c(L,df,pval))
}


###########################################################################
###########################################################################

matern<-function (u, phi, kappa) 
{
    if (is.vector(u)) 
        names(u) <- NULL
    if (is.matrix(u)) 
        dimnames(u) <- list(NULL, NULL)
    uphi <- u/phi
    uphi <- ifelse(u > 0, (((2^(-(kappa - 1)))/ifelse(0, Inf, 
        gamma(kappa))) * (uphi^kappa) * besselK(x = uphi, nu = kappa)), 
        1)
    uphi[u > 600 * phi] <- 0
    return(uphi)
}
########################################################################

dispvar<-function(phi,c0,c1,kap){
set.seed(19650813)
M<-nrow(grid.df)
Tot<-0
for (i in 1:10000){
I<-sample(1:M,1)
J<-sample(1:M,1)
hlag<-distVincentySphere(c(grid.df$Long[I],grid.df$Lat[I]),
c(grid.df$Long[J],grid.df$Lat[J]))/1000 #
Tot<-Tot+c0+c1*(1-matern(hlag, phi, kap))
}
return(Tot/10000)
}

