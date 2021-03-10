###########################################################################
###########################################################################
#
#  R script for exploratory analysis, variogram estimation and validation
#  and kriging of grain nutrient concentrations
#
###########################################################################
###########################################################################


library(xlsx)
library(geosphere)
source("CEPHaStat_2.R")

#  Note  function commands at the bottom of this script must be
#  run before the analyses are attempted
#
#  Note that the script contains options for analysis of the Malawi
#  and the Ethiopia data sets

###################################################################
#
#  Read in the full grain data set
#

read_in.df<-read.xlsx("Ethiopia_Grain.xlsx",sheetIndex=1)
#or
#read_in.df<-read.xlsx("Malawi_Grain.xlsx",sheetIndex=1)

# Select the crop for this analysis
Crop.target<-"Teff"

# Now produce a subset for the target crop

read_in.df<-read_in.df[which(read_in.df$Crop==Crop.target),]

# Now extract element of interest

element<-"Ca"
plotname<-expression("Ca /"~mg~kg^{-1})
plotname2<-expression("Ca /"~log~mg~kg^{-1})

element.col<-which(colnames(read_in.df)==element)

data.df<-read_in.df[,c(2,1,element.col)]

# Extract missing values, if >=77777, 

missing.rows<-which(data.df[,3]>=77777.0)
if(length(missing.rows)>0) {data.df<-data.df[-missing.rows,]}

# Compute summary statistics of the selected variable

summary<-(summa(data.df[,3]))
summary<-rbind(summary,(summa(log(data.df[,3]))))
row.names(summary)<-c("mg/kg,","log mg/kg")

summaplot((data.df[,3]),plotname)
summaplot(log(data.df[,3]),plotname2)
summaplot((outliers(data.df[,3],trim=T)))
summa((outliers(data.df[,3],trim=T)))

# At this point a decision is made whether to analyse the data on their
# original units or to replace them with their natural logarithms.
# In this study the following variables were transformed to logs
#
# Ethiopia: Se and Fe (Teff and Wheat)
# Malawi:  Ca, Se and Fe
#
#  The next two commands are commented out so that they are not 
#  automatically run
#
#    data.df[,3]<-log(data.df[,3])
#    plotname<-plotname2

# The following two rows will remove outliers
# defined as values outwith the outer fences as in Tukey (1977)
# Exploratory Data Analysis, Addison-Wesley.
# This was done for Ethiopia Ca Teff data

#ol<-outliers(data.df[,3])
#data.df<-data.df[-ol,]

#
# Produce a classified post-plot of the data with the bottom
# quartile in blue, the next in green, the next in yellow and
# the top in red.
#

quantiles<-as.character(cut(data.df[,3],quantile(data.df[,3],
c(0,0.25,0.5,0.75,1)),include.lowest=T,
label=c("blue","green","yellow","red")))	


# Classified post-plot
dev.new()
options(scipen=5)
plot(data.df[,1],data.df[,2],pch=16,col=quantiles,asp=1,
xlab="Longitude",ylab="Latitude",cex=0.5)

# Exploratory plots to show E-W or N-S trends

par(mfrow=c(1,2))
par(mar=c(5,5,3,3))
plot(data.df[,1],(data.df[,3]),pch=16,cex=0.8,
xlab="Longitude",ylab=plotname)

plot(data.df[,2],(data.df[,3]),pch=16,cex=0.8,
xlab="Latitude",ylab=plotname)


######################################################################################################################
#
#  Exploratory study of empirical variogram
#
#  Make a variogram cloud
#

N<-nrow(data.df)


NP<-0.5*N*(N-1) 

Long<-data.df$Longitude
Lat<-data.df$Latitude
z<-data.df[,3]

lag<-vector("numeric",NP)
vclo<-vector("numeric",NP)
bear<-vector("numeric",NP)

ico=0

for (i in 1:(N-1)){
for (j in (i+1):N){
ico=ico+1
print(ico/NP)
lag[ico]<-distVincentySphere(c(Long[i],Lat[i]),c(Long[j],Lat[j]))/1000 #
bear[ico]<-finalBearing(c(Long[i],Lat[i]),c(Long[j],Lat[j]), a=6378137, f=1/298.257223563, sphere=TRUE)
zi<-(z[i])
zj<-(z[j])
vclo[ico]<-(zi-zj)
}
}


lagbins<-cut(lag,seq(0,350,35),labels=seq(5,350,35))   # 10-km bins for Ethiopia
#lagbins<-cut(lag,seq(0,100,10),labels=seq(5,100,10))   # 10-km bins for Malawi

lag2<-lag[!is.na(lagbins)]
vclo2<-vclo[!is.na(lagbins)]
lagbins<-factor(lagbins[!is.na(lagbins)])
nlags<-nlevels(lagbins)

####################################################

# Form estimates of the variogram

# Matrices in which to keep, for the lag bins in lagbins:

#  i. semiv - variogram estimates, fifth column is isotropic Matheron, columns
#  1 to 4 are directional (Matheron estimator), columns 6 and 7 are isotropic
#  estimates with Cressie-Hawkins and Dowd estimators respectively.
#
#  ii. npair and varlags - respectively the number of lag pairs per bin, and the mean
#	lag distance.  Columns 1 to 4 are directional (as in semiv) and 5 is isotropic.

semiv<-matrix(nrow=nlags,ncol=7)
npair<-matrix(nrow=nlags,ncol=5)
varlags<-matrix(nrow=nlags,ncol=5)


colnames(semiv)<-c("N-S,Ma","NE-SW,Ma","W-E,Ma","SE-NW,Ma","iso-Ma","iso-CH","iso-Do")
colnames(varlags)<-c("N-S","NE-SW","W-E","SE-NW","iso")
colnames(npair)<-c("N-S","NE-SW","W-E","SE-NW","iso")


varlags[,5]<-as.numeric(by(lag2, lagbins, mean))
npair[,5]<-as.numeric(by(lag2, lagbins, sum))/as.numeric(by(lag2, lagbins, mean))

#
# estimator
#
#  Matheron
semiv[,5]<-0.5*(as.numeric(by((vclo2)^2, lagbins, mean)))

# CH
semiv[,6]<-(as.numeric(by((sqrt(abs(vclo2))), lagbins, mean)))^4
semiv[,6]<-0.5*semiv[,6]/(0.457+0.459/npair[,5]+0.045/npair[,5]^2)

# Do

vclo2d<-c(vclo2,-vclo2)
lagbinsd<-c(lagbins,lagbins)

semiv[,7]<-0.5*(as.numeric(by(vclo2d, lagbinsd, mad)))^2

####################################################

#  Now plot the variograms

maxv<-max(semiv[,5:7])
minv<-min(semiv[,5:7])

plot(varlags[,5],semiv[,5],ylim=c(0,maxv),xlab="Distance /km",ylab="Variance",pch=16)
points(varlags[,5],semiv[,6],pch=1)
points(varlags[,5],semiv[,7],pch=17)

locsv<-minv*c(0.6,0.4,0.2)
points(60,locsv[1],pch=16)
points(60,locsv[2],pch=1)
points(60,locsv[3],pch=17)

text(65,locsv[1],"Matheron",pos=4)
text(65,locsv[2],"Cressie-Hawkins",pos=4)
text(65,locsv[3],"Dowd",pos=4)


####  Fit exponential variogram functions to each set of estimates

# Reasonable initial estimates of the parameters must be entered below

co<-20000
c1<-30000
a<-50

maxl<-max(varlags[,5])*1.1
par(mfrow=c(2,2))

# Matheron's estimator


h<-varlags[,5]
sv<-semiv[,5]
npa<-npair[,5]
kap<-0.5

oo<-optim(c(co,c1,a),wls
,method="L-BFGS-B",
lower=c(0.0,0.0,0),
upper=c(1000000,1000000,5000))
Ahat<-length(h)*log(oo$value)+6  

ooMa<-oo

oon<-optim(c(120),wlsnugg
,method="L-BFGS-B",
lower=c(0.0),
upper=c(1000000))

Ahatn<-length(h)*log(oon$value)+2  

plot(varlags[,5],semiv[,5],ylim=c(0,maxv),xlim=c(0,maxl),xlab="Distance /km",ylab="Variance",pch=16)
#lines(seq(0,200),rep(oon$par[1],201))
lines(seq(0,maxl,0.1),oo$par[1]+oo$par[2]*(1-exp(-seq(0,maxl,0.1)/oo$par[3])))

mtext("a). Matheron",3,adj=0,line=0.5)

# Cressie's and Hawkins's estimator  CH

h<-varlags[,5]
sv<-semiv[,6]
npa<-npair[,5]
kap<-0.5

oo<-optim(c(co,c1,a),wls
,method="L-BFGS-B",
lower=c(0.0,0.0,0),
upper=c(1000000,1000000,5000))
ooCH<-oo
Ahat<-length(h)*log(oo$value)+6  

oon<-optim(c(120),wlsnugg
,method="L-BFGS-B",
lower=c(0.0),
upper=c(1000000))

Ahatn<-length(h)*log(oon$value)+2  


plot(varlags[,5],semiv[,6],ylim=c(0,maxv),xlim=c(0,maxl),xlab="Distance /km",ylab="Variance",pch=16)
lines(seq(0,maxl,0.1),oo$par[1]+oo$par[2]*(1-exp(-seq(0,maxl,0.1)/oo$par[3])))
mtext("b). Cressie-Hawkins",3,adj=0,line=0.5)


# Dowd's estimator Do


h<-varlags[,5]
sv<-semiv[,7]
npa<-npair[,5]
kap<-0.5

oo<-optim(c(co,c1,a),wls
,method="L-BFGS-B",
lower=c(0.0,0.0,0),
upper=c(1000000,1000000,5000))
ooDo<-oo
Ahat<-length(h)*log(oo$value)+6  

oon<-optim(c(120),wlsnugg
,method="L-BFGS-B",
lower=c(0.0),
upper=c(1000000))

Ahatn<-length(h)*log(oon$value)+2  

plot(varlags[,5],semiv[,7],ylim=c(0,maxv),xlim=c(0,maxl),xlab="Distance /km",ylab="Variance",pch=16)
lines(seq(0,maxl,0.1),oo$par[1]+oo$par[2]*(1-exp(-seq(0,maxl,0.1)/oo$par[3])))
mtext("c). Dowd",3,adj=0,line=0.5)

#
#
#  Cross-validation, performed on a subset of the data (all in Ethiopia, 1000 in Malawi)
#
#  An exploratory plot of prediction errors is produced for each variogram
#  Validation outputs are saved in arrays kvopMA, kvopCH and kvopDo for the
#  estimators of Matheron, Cressie & Hawkins and Dowd respectively.
#  

Nv<-N # no need to subset for Ethiopia
#Nv<-500 # Malawi



kvop<-matrix(nrow=Nv,ncol=5)


Nr<-Nv-1

# Make random sample without replacement

sam<-sample(1:N,Nv,replace=F)
#
# Subsample data
Long_s<-Long[sam]
Lat_s<-Lat[sam]
z_s<-(z[sam])

#
# Compute distance matrix for subsample

DM<-matrix(0,nrow=Nv,ncol=Nv)

for (i in 1:(Nv-1)){
for (j in (i+1):Nv){
DM[i,j]<-distVincentySphere(c(Long_s[i],Lat_s[i]),c(Long_s[j],Lat_s[j]))/1000 #WGS84 ellipsoid is used by default
DM[j,i]<-DM[i,j]
}
}

#
#Select parameter set to validate

ooX<-ooMa  # Matheron

# variogram matrix

Gam<-ooX$par[1]+ooX$par[2]*(1-matern(DM,ooX$par[3],0.5))
diag(Gam)<-0

for (i in 1:Nv){

print(i)

A<-Gam[-i,-i]
A<-cbind(A,rep(1,Nr))
A<-rbind(A,c(rep(1,Nr),0))

z_0<-z_s[i]
z_n<-z_s[-i]

b<-c(Gam[i,-i],1)

lam<-solve(A)%*%b
Zhat<-sum(z_n*lam[1:Nr])
err<-z_0-Zhat
kv<-t(b)%*%lam
theta<-(err^2)/kv

kvop[i,]<-c(z_0,Zhat,err,kv,theta)
}

dev.new()
summaplot(kvop[,3],"Kriging error")
kvopMa<-kvop


#Select parameter set to validate

ooX<-ooCH

# variogram matrix

Gam<-ooX$par[1]+ooX$par[2]*(1-matern(DM,ooX$par[3],0.5))
diag(Gam)<-0

for (i in 1:Nv){

print(i)

A<-Gam[-i,-i]
A<-cbind(A,rep(1,Nr))
A<-rbind(A,c(rep(1,Nr),0))

z_0<-z_s[i]
z_n<-z_s[-i]

b<-c(Gam[i,-i],1)

lam<-solve(A)%*%b
Zhat<-sum(z_n*lam[1:Nr])
err<-z_0-Zhat
kv<-t(b)%*%lam
theta<-(err^2)/kv

kvop[i,]<-c(z_0,Zhat,err,kv,theta)
}

dev.new()
summaplot(kvop[,3],"Kriging error")
kvopCH<-kvop


#Select parameter set to validate

ooX<-ooDo

# variogram matrix

Gam<-ooX$par[1]+ooX$par[2]*(1-matern(DM,ooX$par[3],0.5))
diag(Gam)<-0

for (i in 1:Nv){

print(i)

A<-Gam[-i,-i]
A<-cbind(A,rep(1,Nr))
A<-rbind(A,c(rep(1,Nr),0))

z_0<-z_s[i]
z_n<-z_s[-i]

b<-c(Gam[i,-i],1)

lam<-solve(A)%*%b
Zhat<-sum(z_n*lam[1:Nr])
err<-z_0-Zhat
kv<-t(b)%*%lam
theta<-(err^2)/kv

kvop[i,]<-c(z_0,Zhat,err,kv,theta)
}

dev.new()
summaplot(kvop[,3],"Kriging error")
#summa(kvop[,5])
kvopDo<-kvop

colnames(kvopMa)<-c("z_0","Zhat","err","kv","theta")
colnames(kvopCH)<-c("z_0","Zhat","err","kv","theta")
colnames(kvopDo)<-c("z_0","Zhat","err","kv","theta")

#######################################################
#
#  Print out mean and median standardized squared prediction
#  error for each variogram.  If the median for Matheron
#  is in the interval [thetaCLL,thetaCLU] calculated below
#  then it is selected.  Otherwise the variogram for which
#  the median is in the range and closest to 0.455
#

summa(kvopMa[,5])[,c(1,2)]
summa(kvopCH[,5])[,c(1,2)]
summa(kvopDo[,5])[,c(1,2)]

thetaCLU<-0.455+2*sqrt(1/(8*(Nv/2)*(dchisq(0.455,1))^2))
thetaCLL<-0.455-2*sqrt(1/(8*(Nv/2)*(dchisq(0.455,1))^2))

print(c(thetaCLL,thetaCLU))

#Select parameter set for kriging (from ooMa, ooCH, ooDo depending
# on the decision above

ooX<-ooMa

####################################################################
# Writing output prior to kriging.

Models<-rbind(ooMa$par,ooCH$par,ooDo$par)
colnames(Models)<-c("Co","C1","A")
rownames(Models)<-c("Ma","CH","Do")

fname<-paste(element,Crop.target,"fitted_variograms.dat",sep="_")

write.table(Models,fname,quote=F)

fname<-paste(element,Crop.target,"variogram_estimates.dat")

write.table(cbind(varlags,npair,semiv),fname,quote=F,row.names=F)



fname<-paste(element,Crop.target,"XVal_Ma.dat",sep="_")
write.table(kvopMa,fname,quote=F,row.names=F)

fname<-paste(element,Crop.target,"XVal_CH.dat",sep="_")
write.table(kvopCH,fname,quote=F,row.names=F)

fname<-paste(element,Crop.target,"XVal_Do.dat",sep="_")
write.table(kvopDo,fname,quote=F,row.names=F)

fname<-paste(element,Crop.target,"XVal_sample.dat",sep="_")
write.table(sam,fname,quote=F,row.names=F)

fname<-paste(element,Crop.target,"_summary.dat",sep="_")
write.table(summary,fname,quote=F)

#####################################################################
#
#  Kriging
#
targets.df<-read.csv("Ethiopia_map.csv",header=T) #Option for Ethiopia
#targets.df<-read.csv("Malawi_targets.csv",header=T) #Option for Malawi

urdata.df<-data.frame(cbind(Lat,Long,z))

# Make a distance matrix among all data points 

Long<-urdata.df$Long
Lat<-urdata.df$Lat

DM<-matrix(0,nrow=N,ncol=N)

for (i in 1:(N-1)){
for (j in (i+1):N){

DM[i,j]<-distVincentySphere(c(Long[i],Lat[i]),c(Long[j],Lat[j]))/1000 #WGS84 ellipsoid is used by default
DM[j,i]<-DM[i,j]
}
}

# variogram matrix

Gam<-ooX$par[1]+ooX$par[2]*(1-matern(DM,ooX$par[3],0.5))

A<-Gam
A<-cbind(A,rep(1,N))
A<-rbind(A,c(rep(1,N),0))
Ainv<-solve(A)


Nt<-nrow(targets.df)

krop<-matrix(nrow=Nt,ncol=4)

Allpoints<-cbind(Lat,Long)


for (it in 1:Nt){
print(c(it,Nt))

# Extract Lat and Long of target point 

Lat_t<-targets.df$Lat[it]
Long_t<-targets.df$Long[it]

Bd<-matrix(nrow=(N+1),ncol=1)

for (j in 1:N){
Bd[j]<-distVincentySphere(c(Long[j],Lat[j]),c(Long_t,Lat_t))/1000
}

#maxdis<-max(Bd[1:Ne])
#print(c(it,it/Nt,maxdis))

Bd[1:N]<-ooX$par[1]+ooX$par[2]*(1-matern(Bd[1:N],ooX$par[3],0.5))
Bd[N+1]<-1
lam<-Ainv%*%Bd
Zhat<-sum(z*lam[1:N])
kv<-t(Bd)%*%lam
lagr<-lam[N+1]
krop[it,]<-c(Long_t,Lat_t,Zhat,kv)
}


colnames(krop)<-c("Longitude","Latitude","Zhat","kv")

fname<-paste(element,"_OK_",Crop.target,".csv",sep="")

write.csv(krop,fname,row.names=F)  



###################################################################################################
#
# FUNCTIONS
#
###################################################################################################

############################################################################
#
#  Matern correlation function
#
#

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


############################

wls<-function(theta){

c0<-theta[1]
c1<-theta[2]
phi<-theta[3]


nh<-length(h)

wss<-0

for (i in 1:nh){
svm<-c0+c1*(1-matern(h[i],phi,kap))
wss<-wss+npa[i]*((sv[i]-svm)^2)
}
wss<-wss/nh
return(wss)
}


wlsnugg<-function(theta){

c0<-theta[1]

nh<-length(h)

wss<-0

for (i in 1:nh){
svm<-c0
wss<-wss+npa[i]*((sv[i]-svm)^2)
}
wss<-wss/nh
return(wss)
}






