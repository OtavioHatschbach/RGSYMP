library(survival)

kmquantiles<-function(time, status, x, p=0.75, points=101){
	
	b<-bw.SJ(x)
	
	xs<-seq(min(x),max(x),length=points)
	rval <- list(x=xs, y=sapply(xs,function(xi) kmq(time,status,x,p,xi,b)))
	if(length(p)>1)
		rownames(rval$y)<-p
	rval
}

kmq<-function(time,status,x,p, x0,b){
	w<-dnorm(x,x0,s=b)
	S<-survfit(Surv(time,status)~1,weights=w)
    sapply(p, function(p.){
    		i<-suppressWarnings(min(which(S$surv<=p.)))
    		S$time[i]
    })
}

kmgrid<-function(time, status, x, nx=61,nt=51){
	xgrid<-seq(min(x),max(x),length=nx)
	tgrid<-seq(min(time),max(time),length=nt)
	
	b<-bw.SJ(x)
	Sgrid<-matrix(ncol=nt,nrow=nx)
	for(i in 1:nx){
		w<-dnorm(x,xgrid[i],s=b)
		S<-survfit(Surv(time,status)~1,weights=w)
		Sgrid[i,]<-approx(S$time, S$surv, tgrid, method="constant",yleft=1,rule=2,f=1)$y
	}
	list(x=xgrid, y=tgrid, z=Sgrid)
}

plotSx<-function(time, status, x,x0,...){
		b<-bw.SJ(x)
		w<-dnorm(x,x0,s=b)
		S<-survfit(Surv(time,status)~1,weights=w)
		plot(S,...)
}

plotROCt<-function(time,status,x,t0,nx=61,type="s",xlab="1-Spec",ylab="Sens",...){
	b<-bw.SJ(x)

	xgrid<-seq(min(x),max(x),length=nx)
	Sw<-function(z){
		w<-dnorm(x,z,s=b)
		S<-survfit(Surv(time,status)~1,weights=w)
		min(S$surv[S$time<=t0])
	}
   Sgrid<-sapply(xgrid,Sw)
   n<-length(time)
   pgrid<-diff(c(0,sapply(xgrid, function(z) sum(x<=z))))
   
   Scond<-cumsum(rev(Sgrid*pgrid))/n
   mFx<-cumsum(rev(pgrid))/n
   sens<-(mFx-Scond)/(1-Scond[nx])
   mspec<-Scond/Scond[nx]
   plot(mspec, sens, type=type,xlab=xlab,ylab=ylab,...)
}

hatdotplot<-function(time,status, x,t0,mean=TRUE){
	b<-bw.SJ(x)

	Sw<-function(z){
		w<-dnorm(x,z,s=b)
		S<-survfit(Surv(time,status)~1,weights=w)
		min(S$surv[S$time<=t0])
	}
	hatS<-sapply(x, Sw)
	xdead<-x
	xalive<-x
	pdead<-1-hatS
	palive<-hatS
	plot(c(xdead,xalive), jitter(rep(1:0,each=length(x))), 
		pch=19,cex=1.5*sqrt(c(pdead,palive)),
		col="#00000040",ylab="P(Dies)",xlab=deparse(substitute(x)),
		yaxt="n")
	axis(2, at=c(0,1))
	if(mean){
		segments(y0=0.5,x0=weighted.mean(xdead,pdead),y1=1.5,col="red",lwd=2)
		segments(y=-0.5,x0=weighted.mean(xalive,palive),y1=0.5,col="red",lwd=2)
	}
}




triplot<-function(time,status,x){
	quartz(height=5,width=5)
  p<-NULL
	repeat({
		dev.set(2)
		plot(time~x,pch=ifelse(status==0,1,19),
	 		col=c("grey","orange","darkred")[status+1], 
	 		xlab="log(prothrombin time)",ylab="Time (years)")
		if (!is.null(p)){
		  abline(h=p$y,lty=1,lwd=2,col="forestgreen")
		  abline(v=p$x,lty=1,lwd=2,col="forestgreen")
		  
		}
		p<-locator(1)
		if(is.null(p)) break
		
		abline(h=p$y,lty=1,lwd=2,col="green")
		abline(v=p$x,lty=1,lwd=2,col="green")

		dev.set(3)
		par(mfrow=c(2,1),pty="m",mar=c(4.1,2.1,.1,.1))
	  	plotSx(time,status>0,x,p$x)
	  	par(pty="s")
  		plotROCt(time,status>0,x,p$y)
  		abline(0,1,lty=3,col="grey")
	})
dev.off(3)
}

triplot2<-function(time,status,x){
	quartz(height=5,width=5)
  p<-NULL
	repeat({
		dev.set(2)
		plot(time~x,pch=ifelse(status==0,1,19),
	 		col=c("grey","orange","darkred")[status+1], 
	 		xlab=deparse(substitute(x)),ylab="Time (years)")
		if (!is.null(p)){
		  abline(h=p$y,lty=1,lwd=2,col="forestgreen")
		  abline(v=p$x,lty=1,lwd=2,col="forestgreen")
		}		
		p<-locator(1)
		if(is.null(p)) break
		
		abline(h=p$y,lty=1,lwd=2,col="green")
		abline(v=p$x,lty=1,lwd=2,col="green")

		dev.set(3)
		par(mfrow=c(2,1),pty="m",mar=c(4.1,2.1,.1,.1))
	  	plotSx(time,status>0,x,p$x)
  		hatdotplot(time,status>0,x,p$y)
	})
dev.off(3)
}

recolor<-function(S){
	nrz<-51
	ncz<-61
	z<-t(Sxt$z)
	jet.colors <- colorRampPalette( c("blue", "green") )
	# Generate the desired number of colors from this palette
	nbcol <- 100
	color <- jet.colors(nbcol)
	# Compute the z-value at the facet centres
	zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
	# Recode facet z-values into color indices
	facetcol <- cut(zfacet, nbcol)
	color[facetcol]
}