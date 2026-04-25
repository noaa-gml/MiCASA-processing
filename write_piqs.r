#!/usr/bin/env Rscript

ct.setup()

din <- load.ncdf("monthly_1x1/MiCASA_v1_flux_x360_y180_monthly.nc")

gpp <- -2*din$NPP
rtot <- din$Rh + din$NPP

nee1 <- din$Rh-din$NPP
nee2 <- gpp+rtot

nmon <- length(din$time)
plt.start <- as.POSIXlt(din$time[1])
y0 <- plt.start$year+1900
m0 <- plt.start$mon+1
x.time <- as.numeric(seq(ISOdatetime(y0,m0,1,0,0,0,tz="UTC"),
                         by="1 month",length.out=nmon+1))



plot.tser.cubic <- function(x,ybar,fit,main='') {
  par(las=1)
  n <- length(ybar)
  plot(NA,NA,xlim=range(x),ylim=c(0.8,1.2)*range(fit$y),xlab='x',ylab='y',
       xaxt='n',main=main)
  abline(v=0)
  axis.POSIXct(side=1,x=x,at=seq(ISOdatetime(2000,1,1,0,0,0,tz="UTC"),
                            ISOdatetime(2016,1,1,0,0,0,tz="UTC"),
                            by="3 months"))
  for( i in 1:n) {
    from=x[i]
    to=x[i+1]
    curve(fit$a[i]*(x-from)^3+fit$b[i]*(x-from)^2+fit$c[i]*(x-from)+fit$d[i],from=from,to=to,add=TRUE)
    abline(v=to,lty=3)
    lines(x=c(from,to),y=rep(ybar[i],2),col="blue")
    points(from,fit$y[i])
  }
  points(to,fit$y[n+1])
}


plot.tser.quadratic <- function(x,ybar,fit,main='') {
  par(las=1)
  n <- length(ybar)
  plot(NA,NA,xlim=range(x),ylim=c(0.8,1.2)*range(fit$y),xlab='x',ylab='y',
       xaxt='n',main=main)
  abline(v=0)
  axis.POSIXct(side=1,x=x,at=seq(ISOdatetime(2000,1,1,0,0,0,tz="UTC"),
                            ISOdatetime(2016,1,1,0,0,0,tz="UTC"),
                            by="3 months"))
  par(xpd=NA)
  for( i in 1:n) {
    from=x[i]
    to=x[i+1]
    curve(fit$a[i]*(x-from)^2+fit$b[i]*(x-from)+fit$c[i],from=from,to=to,add=TRUE)
    abline(v=to,lty=3)
    lines(x=c(from,to),y=rep(ybar[i],2),col="blue")
    points(from,fit$y[i])
  }
  points(to,fit$y[n+1])
}


plot.tser.linear <- function(x,ybar,fit,main='') {
  par(las=1)
  n <- length(ybar)
  plot(NA,NA,xlim=range(x),ylim=c(0.8,1.2)*range(fit$y),xlab='x',ylab='y',
       xaxt='n',
       main=main)
  abline(v=0)
  axis.POSIXct(side=1,x=x,at=seq(ISOdatetime(2000,1,1,0,0,0,tz="UTC"),
                            ISOdatetime(2016,1,1,0,0,0,tz="UTC"),
                            by="3 months"))
  for( i in 1:n) {
    from=x[i]
    to=x[i+1]
    lines(x[i:(i+1)],
          fit$y[i:(i+1)])
    abline(v=to,lty=3)
    lines(x=c(from,to),y=rep(ybar[i],2),col="blue")
    points(from,fit$y[i])
  }
  points(to,fit$y[n+1])
}

pb <- progress.bar.start(360*180,360*180)

piqsfit.gpp <- list()
piqsfit.resp <- list()

piqsfit.gpp$a <- array(NA,dim=c(360,180,nmon))
piqsfit.gpp$b <- array(NA,dim=c(360,180,nmon))
piqsfit.gpp$c <- array(NA,dim=c(360,180,nmon))
#piqsfit.gpp$d <- array(NA,dim=c(360,180,nmon))

piqsfit.resp$a <- array(NA,dim=c(360,180,nmon))
piqsfit.resp$b <- array(NA,dim=c(360,180,nmon))
piqsfit.resp$c <- array(NA,dim=c(360,180,nmon))
#piqsfit.resp$d <- array(NA,dim=c(360,180,nmon))

ipb <- 0
for (i in 1:360) {
  for (j in 1:180) {
    ipb <- ipb+1
    
    if(all(nee1[i,j,] == 0)) {
      next
    }

    fit.gpp <- piqs(x.time,gpp[i,j,])
    fit.rtot <- piqs(x.time,rtot[i,j,])
    #    fit.gpp <- pils.2(x.time,gpp[i,j,])
    #    fit.rtot <- pils.2(x.time,rtot[i,j,])
    #    fit.gpp <- pics(x.time,gpp[i,j,])
    #    fit.rtot <- pics(x.time,rtot[i,j,])
    if(FALSE) {
      pdf(file=sprintf("piqs_i%d_j%d.pdf",i,j),width=15,height=8)
      layout(matrix(1:2,2,1))
      plot.tser.quadratic(x=x.time,ybar=gpp[i,j,],fit=fit.gpp,main=sprintf("GPP at ij=(%d,%d), lon %.1f, lat %.1f",i,j,din$longitude[i],din$latitude[j]))
      
      plot.tser.quadratic(x=x.time,ybar=rtot[i,j,],fit=fit.rtot,main=sprintf("Rtot at ij=(%d,%d), lon %.1f, lat %.1f",i,j,din$longitude[i],din$latitude[j]))
      #    plot.tser.cubic(x=x.time,ybar=gpp[i,j,],fit=fit.gpp,main=sprintf("GPP at ij=(%d,%d), lon %.1f, lat %.1f",i,j,din$lon[i],din$lat[j]))
      
      #    plot.tser.cubic(x=x.time,ybar=rtot[i,j,],fit=fit.rtot,main=sprintf("Rtot at ij=(%d,%d), lon %.1f, lat %.1f",i,j,din$lon[i],din$lat[j]))
      dev.off()
    }

    piqsfit.gpp$a[i,j,] <- fit.gpp$a
    piqsfit.gpp$b[i,j,] <- fit.gpp$b
    piqsfit.gpp$c[i,j,] <- fit.gpp$c
#    piqsfit.gpp$d[i,j,] <- fit.gpp$d
    
    piqsfit.resp$a[i,j,] <- fit.rtot$a
    piqsfit.resp$b[i,j,] <- fit.rtot$b
    piqsfit.resp$c[i,j,] <- fit.rtot$c
#    piqsfit.resp$d[i,j,] <- fit.rtot$d
    
    pb <- progress.bar.print(pb,ipb)
    
  }
}

progress.bar.end(pb)

piqsfit.time <- x.time[1:(length(x.time)-1)]
save(file='fit.piqs.rda',piqsfit.gpp,piqsfit.resp,piqsfit.time)
