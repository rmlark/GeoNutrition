#############################################################################
#
#  CEPHaStat:  R functions produced for the CEPHaS GCRF Project, Working 
#	Group 4 for use in CEPHaS and related research on conservation 
#	agriculture
#
vers<-"v1.1 Correction to order of summa output"
#

header<-"\n\nCEPHaStat is a collection of R functions produced by the UKRI GCRF
CEPHaS project.  The functions relate primarily to statistical analysis of
soil and related data from agricultural and environmental surveys and 
experiments.  The CEPHaStat file can be shared as it stands, and is made 
available to all for use, without warranty or liability.\n\n"
#
fnclist<-"Current list of functions\n
skew kurt ocskew
summa summaplot outliers
okgridvar, sv, dis
vanG, vanGdraw, campB, campDdraw, plotgroups\n \n"
#############################################################################

cat(paste("CEPHaStat version ",vers))
cat(header)
cat(fnclist)
#############################################################################

skew<-function(x){

# compute coefficient of skewness of data in x

x<-na.drop(x)
n<-length(x)
xd<-x-mean(x)
mu3<-sum(xd^3)/(n-1)
mu<-sqrt(sum(xd^2)/(n-1))
sk<-mu3/(mu^3)
return(sk)
}

#############################################################################

kurt<-function(x){

# compute coefficient of kurtosis of data in x

x<-na.drop(x)
n<-length(x)
xd<-x-mean(x)
mu4<-sum(xd^4)/(n-1)
mu<-sqrt(sum(xd^2)/(n-1))
sk<-(mu4/(mu^4))-3
return(sk)
}

################################################################################

ocskew<-function(x){

# compute the octile skewness of values in x

x<-na.drop(x)
Ocs<-quantile(x,c(1/8,0.5,7/8))
os<-((Ocs[3]-Ocs[2])-(Ocs[2]-Ocs[1]))/(Ocs[3]-Ocs[1])
return(os)
}

###########################################################################

summa<-function(x,sigf){

# compute summary statistics of values in x

if(missing(sigf)){rosig<-F}else{rosig<-T}


x<-na.drop(x)
Q1<-quantile(x,prob=0.25)
Q3<-quantile(x,prob=0.75)
Q2<-quantile(x,prob=0.5)
hspread<-Q3-Q1
Fu<-Q3+3*hspread
Fl<-Q1-3*hspread

# ols,oll: values below and above outer fences
# posols,posoll: values below and above inner fences 
# (so ols and posols overlap, as do oll and posoll
#

ols<-which(x<Fl)
oll<-which(x>Fu)
posols<-which(x<(Q1-1.5*hspread))
if(length(posols)==0){
lw<-min(x)}else{
lw<-min(x[-posols])}
posoll<-which(x>(Q3+1.5*hspread))
if(length(posoll)==0){
uw<-max(x)}else{
uw<-max(x[-posoll])}

ol<-c(ols,oll) # combined outlier set

nol<-length(ol)

outp<-matrix(c(mean(x),Q2,Q1,Q3,var(x),sqrt(var(x)),skew(x),ocskew(x),kurt(x),nol),1,10)

if(rosig=="TRUE"){outp<-signif(outp,sigf)}

colnames(outp)<-c("Mean","Median",
"Quartile.1", "Quartile.3","Variance","SD","Skewness",
"Octile skewness","Kurtosis",
"No. outliers")


return(outp)

}

####################################################################

summaplot<-function(x,varname){

# plot a histogram with boxplot and QQ plot of data in x indicating
# any probable outliers by Tukey's criterion

x<-na.drop(x)
if(missing(varname)){varname<-"x"}

Q1<-quantile(x,prob=0.25)
Q3<-quantile(x,prob=0.75)
Q2<-quantile(x,prob=0.5)
hspread<-Q3-Q1
Fu<-Q3+3*hspread
Fl<-Q1-3*hspread


# ols,oll: values below and above outer fences
# posols,posoll: values below and above inner fences 
# (so ols and posols overlap, as do oll and posoll
#
ols<-which(x<Fl)
oll<-which(x>Fu)
posols<-which(x<(Q1-1.5*hspread))
if(length(posols)==0){
lw<-min(x)}else{
lw<-min(x[-posols])}
posoll<-which(x>(Q3+1.5*hspread))
if(length(posoll)==0){
uw<-max(x)}else{
uw<-max(x[-posoll])}

ol<-c(ols,oll) # combined outlier set
par(mfrow=c(1,2))
ymax<-max((hist(x,plot=F))$counts)
hist(x,main="",col="AliceBlue", xlab=varname,ylim=c(0,(ymax*1.25)))

boxmin<-ymax*1.1
boxmax<-ymax*1.2
boxmid<-ymax*1.15

lines(c(Q1,Q3),c(boxmin,boxmin))
lines(c(Q1,Q3),c(boxmax,boxmax))
lines(c(Q1,Q1),c(boxmin,boxmax))
lines(c(Q3,Q3),c(boxmin,boxmax))
lines(c(Q1,lw),c(boxmid,boxmid))
lines(c(Q3,uw),c(boxmid,boxmid))
lines(c(Q2,Q2),c(boxmin,boxmax),lwd=2)

lines(c(Fu,Fu),c(10,boxmid),lty=5,col="red")
lines(c(Fl,Fl),c(10,boxmid),lty=5,col="red")

qqn<-qqnorm(x,main="",pch=16)
qqline(x)
points(qqn$x[ol],qqn$y[ol],pch=16,col="red")

}
####################################################################

histplot<-function(x,varname){

# plot a histogram with boxplot in x indicating
# any probable outliers by Tukey's criterion

x<-na.drop(x)
if(missing(varname)){varname<-"x"}

Q1<-quantile(x,prob=0.25)
Q3<-quantile(x,prob=0.75)
Q2<-quantile(x,prob=0.5)
hspread<-Q3-Q1
Fu<-Q3+3*hspread
Fl<-Q1-3*hspread


# ols,oll: values below and above outer fences
# posols,posoll: values below and above inner fences 
# (so ols and posols overlap, as do oll and posoll
#
ols<-which(x<Fl)
oll<-which(x>Fu)
posols<-which(x<(Q1-1.5*hspread))
if(length(posols)==0){
lw<-min(x)}else{
lw<-min(x[-posols])}
posoll<-which(x>(Q3+1.5*hspread))
if(length(posoll)==0){
uw<-max(x)}else{
uw<-max(x[-posoll])}

ol<-c(ols,oll) # combined outlier set

ymax<-max((hist(x,plot=F))$counts)
hist(x,main="",col="AliceBlue", xlab=varname,ylim=c(0,(ymax*1.25)))

boxmin<-ymax*1.1
boxmax<-ymax*1.2
boxmid<-ymax*1.15

lines(c(Q1,Q3),c(boxmin,boxmin))
lines(c(Q1,Q3),c(boxmax,boxmax))
lines(c(Q1,Q1),c(boxmin,boxmax))
lines(c(Q3,Q3),c(boxmin,boxmax))
lines(c(Q1,lw),c(boxmid,boxmid))
lines(c(Q3,uw),c(boxmid,boxmid))
lines(c(Q2,Q2),c(boxmin,boxmax),lwd=2)

lines(c(Fu,Fu),c(10,boxmid),lty=5,col="red")
lines(c(Fl,Fl),c(10,boxmid),lty=5,col="red")

}
###########################################################################

outliers<-function(x,trim){

# compute summary statistics of values in x

if(missing(trim)){trim<-F}

x<-na.drop(x)
Q1<-quantile(x,prob=0.25)
Q3<-quantile(x,prob=0.75)
Q2<-quantile(x,prob=0.5)
hspread<-Q3-Q1
Fu<-Q3+3*hspread
Fl<-Q1-3*hspread

# ols,oll: values below and above outer fences
# posols,posoll: values below and above inner fences 
# (so ols and posols overlap, as do oll and posoll
#

ols<-which(x<Fl)
oll<-which(x>Fu)
posols<-which(x<(Q1-1.5*hspread))
if(length(posols)==0){
lw<-min(x)}else{
lw<-min(x[-posols])}
posoll<-which(x>(Q3+1.5*hspread))
if(length(posoll)==0){
uw<-max(x)}else{
uw<-max(x[-posoll])}

ol<-c(ols,oll) # combined outlier set

nol<-length(ol)

print(paste(nol," value(s) lie outwith the outer fences (Tukey, 1977)"),quote=F)

if(nol!=0){
print(paste("Outlying values are:"),quote=F)
print(x[ol])
print(paste("Indices of outlying values are"),quote=F)

if(trim==F){
return(ol)}else{
print(ol)
return(x[-ol])
}
}
}


###############################################################
###############################################################

okgridvar<-function(space,modtype,c0,c1,a1){

x0<-c(0,0)
x<-rep(seq(-2.5,2.5,1),6)*space
y<-rep(seq(-2.5,2.5,1),rep(6,6))*space
xg<-cbind(x,y)

a<-matrix(1,37,37)
a[37,37]<-0
b<-rep(1,37)

for (i in 1:36){
lag<-ddis(i,xg)
b[i]<-sv(lag,modtype,c0,c1,a1)
for (j in i:36){
if(i==j){a[i,j]<-0.0
}else{
lag<-dis(i,j,xg)
a[i,j]<-sv(lag,modtype,c0,c1,a1)
a[j,i]<-a[i,j]}
}
}

ai<-solve(a)
lam<-ai%*%b
kv<-t(b)%*%lam
return(kv)
}

###############################################################
###############################################################
#
#  Functions below not yet documented
#
###############################################################
###############################################################


sv<-function(lag,modtype,c0,c1,a1){

hovera<-lag/a1

#modtype<-"Sph" or "Exp" 

if(modtype=="Exp"){semiv<-c0+c1*(1-exp(-hovera))
} else { 
if(lag>a1){semiv<-c0+c1}
if(lag<=a1){sph<-1.5*hovera-0.5*hovera^3
semiv<-c0+c1*sph}
}
return(semiv)
}

###############################################################
###############################################################

dis<-function(i,j,xg){

dx<-xg[i,1]-xg[j,1]
dy<-xg[i,2]-xg[j,2]
return(sqrt((dx*dx)+(dy*dy)))
}

ddis<-function(i,xg){

dx<-xg[i,1]
dy<-xg[i,2]
return(sqrt((dx*dx)+(dy*dy)))
}

ocskew<-function(x){

percs<-quantile(x,c(1/8,1/2,7/8))
o7<-percs[3]
o4<-percs[2]
o1<-percs[1]

os<-((o7-o4)-(o4-o1))/(o7-o1)
return(as.numeric(os))
}


#####################################################################
na.drop<-function(xin){
noNA<-as.numeric(length(which(is.na(xin)==T)))
if(noNA>0){
x<-as.numeric(na.omit(xin))
print(paste(noNA," missing value(s) removed"),quote=F)
}else{
x<-xin
}
return(x)
}

######################################################################


vanG<-function(ps1,id,xidep){
h<-xidep[,1]
thr<-ps1[id,1]
ths<-ps1[id,2]
alp<-ps1[id,3]
nscal<-ps1[id,4]
mscal<-ps1[id,5]
den<-(1+((abs(alp*h))^nscal))^mscal
return(thr+((ths-thr)/den))
}

vanGdraw<-function(psi,h){
thr<-psi[1]
ths<-psi[2]
alp<-psi[3]
nscal<-psi[4]
mscal<-psi[5]
den<-(1+((abs(alp*h))^nscal))^mscal
return(thr+((ths-thr)/den))
}

campB<-function(ps1,id,xidep){
h<-xidep[,1]
ths<-ps1[id,1]
alp<-ps1[id,2]
nscal<-ps1[id,3]
return(ths*((alp*h)^(-1/nscal)))
}


campBdraw<-function(psi,h){

ths<-psi[1]
alp<-psi[2]
nscal<-psi[3]

return(ths*((alp*h)^(-1/nscal)))
}



plotgroups<-function(x,y,icl,data.df,xlab,ylab,xylog){
# function to plot variables x and y in dataframe data.df
# groups with different variables (in icl) are plotted with contrasting
# symbols.  xlab and ylab are optional labels.  Set xylog to x or
# y or xy for log scales, n for neither

options(scipen=999)
if(missing(xlab)){xlab<-"x"}
if(missing(ylab)){ylab<-"y"}
if(missing(xylog)){xylog<-"n"}

cnames<-colnames(data.df)
indl<-which(cnames==icl)
indx<-which(cnames==x)
indy<-which(cnames==y)

levels(data.df[,indx])<-seq(1,nlevels(data.df[,indx]))


if(xylog=="n"){
plot(data.df[,indx],data.df[,indy],xlab=xlab,ylab=ylab,pch=16,col=data.df[,indl])
}else{
plot(data.df[,indx],data.df[,indy],log=xylog,xlab=xlab,ylab=ylab,pch=16,col=data.df[,indl])
}

}


mean_censor<-function(y,cen,log.t){
if(missing(log.t)){log.t=F}

mean.guess<-mean(y)
sd.guess<-sd(y)

ncen<-length(which(is.na(censor(y,cen))))
oo<-optim(c(mean.guess,sd.guess),nllcen,cen=cen,y=y)
oo2<-optim(oo$par,nllcen,cen=cen,y=y,hessian=T)

se<-sqrt((solve(oo2$hessian))[1,1])

mean.est<-oo2$par[1]
sd.est<-oo2$par[2]

# approximate df by n-1

nobs<-length(y)
tval<-qt(0.975,(nobs-1))
clu<-oo2$par[1]+tval*se
cll<-oo2$par[1]-tval*se

mean.naive<-mean(censor(y,cen),na.rm=T)

op1<-matrix(c((-1*oo2$value),ncen,mean.est,sd.est,cll,clu,mean.naive),1,7)
colnames(op1)<-c("log_likelihood","Number_censored","Mean","SD","Lower_95%","Upper_95%","Mean_>_DL_only")



if(log.t==T){

median.unbiased<-exp(mean.est)
cll.bt<-exp(cll)
clu.bt<-exp(clu)


op2<-matrix(c(median.unbiased,cll.bt,clu.bt),1,3)
colnames(op2)<-c("Back_median","Back_Lower_95%","Back_Upper_95%")
op1<-cbind(op1,op2)
}

return(op1)
}

#
censor<-function(x,cen){
censored<-which(x<cen)
x[censored]<-NA
return(x)
}

nllcen<-function(theta,cen,y){
#
#
#
mu<-theta[1]
sig<-theta[2]

y<-censor(y,cen)
n_c<-length(which(is.na(y)))
Phi_c<-log(pnorm(cen,mu,sig))
z<-y[-which(is.na(y))]

dens<-dnorm(z,mu,sig,log=T)

return(-n_c*Phi_c-sum(dens))
}



