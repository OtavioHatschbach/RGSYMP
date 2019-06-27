
source("rg-code.R")
data(pbc)
pbc$timeyrs<-pbc$time/365.25
pbc<-subset(pbc, !is.na(trt))
plot(timeyrs~log(protime),xlab="log(prothrombin time)",ylab="Time (years)",data=pbc)

plot(timeyrs~log(protime),pch=ifelse(status==0,1,19), col=c("grey","orange","darkred")[status+1],xlab="log(prothrombin time)",ylab="Time (years)",data=pbc)

qq<-with(pbc, kmquantiles(timeyrs, status>0, log(protime),p=c(0.95,0.9,0.75,0.5)))
matlines(qq$x, t(qq$y),lty=1,col="black")



Sxt<-with(pbc, kmgrid(timeyrs, status>0, log(protime)))
with(Sxt, image(x,y,z,xlab="log(protime)",ylab="Time(years)",main="Survival"))

pbc$protime[pbc$id==107]<-10.7
Sxt<-with(pbc, kmgrid(timeyrs, status>0, log(protime)))
with(Sxt, image(x,y,z,xlab="log(protime)",ylab="Time(years)",main="Survival"))



with(Sxt, filled.contour(x,y,z,xlab="log(protime)",ylab="Time(years)",main="Survival"))

opar<-par(mar=c(0,0,0,0))
with(Sxt, persp(y,x,t(z),ylab="log(protime)",xlab="Time(years)",zlab="Survival",theta=110,phi=30,shade=0.6,expand=0.5))

colours<-recolor(Sxt)

with(Sxt, persp(y,x,t(z), col = colours,ylab="log(protime)",xlab="Time(years)",zlab="Survival",theta=110,phi=30,shade=0.6,expand=0.5,border=NA))

par(opar)

opar<-par(mar=c(4.1,2.1,0.1,0.1))
with(pbc, triplot(timeyrs, status, log(protime)))

par(mar=c(4.1,2.1,0.1,0.1),mfrow=c(1,1))
with(pbc, triplot2(timeyrs, status, log(protime)))

with(pbc, triplot2(timeyrs, status, log(bili)))


