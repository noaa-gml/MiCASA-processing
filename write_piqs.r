#!/usr/bin/env Rscript

ct.setup()
source(file.path(Sys.getenv("WORK_DIR", getwd()), "config.r"))
cfg <- micasa.config()

din <- load.ncdf(micasa.out.monthly.cat(cfg))

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

# ---------------------------------------------------------------------------
# Edge padding (proposal #1 in README.ash)
#
# Default (no env vars set): script behaves exactly as before. Set
# MICASA_PIQS_PAD_RIGHT (and optionally MICASA_PIQS_PAD_LEFT) to a small
# positive integer to pad the time series with that many synthetic months
# before fitting. The pad values are the per-cell climatology of the same
# calendar month, drawn from the unpadded data. Pad coefficients are stripped
# before saving, so piqsfit.gpp/resp arrays and piqsfit.time keep their
# original shape.
#
# Recommended starting point in production: PAD_RIGHT=2, PAD_LEFT=0.
# ---------------------------------------------------------------------------
pad.left  <- as.integer(Sys.getenv("MICASA_PIQS_PAD_LEFT",  unset="0"))
pad.right <- as.integer(Sys.getenv("MICASA_PIQS_PAD_RIGHT", unset="0"))
if(is.na(pad.left)  || pad.left  < 0) pad.left  <- 0
if(is.na(pad.right) || pad.right < 0) pad.right <- 0
cat(sprintf("PIQS edge padding: left=%d, right=%d (set MICASA_PIQS_PAD_{LEFT,RIGHT} to override)\n",
            pad.left, pad.right))

pad.start.lt <- as.POSIXlt(ISOdatetime(y0, m0, 1, 0, 0, 0, tz="UTC"))
pad.start.lt$mon <- pad.start.lt$mon - pad.left
x.time.ext <- as.numeric(seq(as.POSIXct(pad.start.lt, tz="UTC"),
                             by="1 month",
                             length.out=nmon + pad.left + pad.right + 1))

data.month <- as.POSIXlt(as.POSIXct(x.time[1:nmon], origin="1970-01-01", tz="UTC"))$mon + 1
seg.month  <- as.POSIXlt(as.POSIXct(x.time.ext[1:(nmon + pad.left + pad.right)],
                                    origin="1970-01-01", tz="UTC"))$mon + 1
pad.idx <- integer(0)
if(pad.left  > 0) pad.idx <- c(pad.idx, seq_len(pad.left))
if(pad.right > 0) pad.idx <- c(pad.idx, seq(nmon + pad.left + 1, nmon + pad.left + pad.right))



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

    # When pad.left == pad.right == 0 this is a no-op and the call to piqs()
    # is identical to the original code path.
    if(pad.left > 0 || pad.right > 0) {
      gpp.cell.ext  <- numeric(nmon + pad.left + pad.right)
      rtot.cell.ext <- numeric(nmon + pad.left + pad.right)
      keep <- (pad.left + 1):(pad.left + nmon)
      gpp.cell.ext[keep]  <- gpp[i,j,]
      rtot.cell.ext[keep] <- rtot[i,j,]
      for(p in pad.idx) {
        cm <- seg.month[p]
        same.cm <- which(data.month == cm)
        gpp.cell.ext[p]  <- mean(gpp[i,j,same.cm],  na.rm=TRUE)
        rtot.cell.ext[p] <- mean(rtot[i,j,same.cm], na.rm=TRUE)
      }
    } else {
      gpp.cell.ext  <- gpp[i,j,]
      rtot.cell.ext <- rtot[i,j,]
      keep <- 1:nmon
    }

    fit.gpp <- piqs(x.time.ext, gpp.cell.ext)
    fit.rtot <- piqs(x.time.ext, rtot.cell.ext)
    #    fit.gpp <- pils.2(x.time.ext, gpp.cell.ext)
    #    fit.rtot <- pils.2(x.time.ext, rtot.cell.ext)
    #    fit.gpp <- pics(x.time.ext, gpp.cell.ext)
    #    fit.rtot <- pics(x.time.ext, rtot.cell.ext)
    if(FALSE) {
      pdf(file=sprintf("piqs_i%d_j%d.pdf",i,j),width=15,height=8)
      layout(matrix(1:2,2,1))
      plot.tser.quadratic(x=x.time.ext,ybar=gpp.cell.ext,fit=fit.gpp,main=sprintf("GPP at ij=(%d,%d), lon %.1f, lat %.1f",i,j,din$longitude[i],din$latitude[j]))

      plot.tser.quadratic(x=x.time.ext,ybar=rtot.cell.ext,fit=fit.rtot,main=sprintf("Rtot at ij=(%d,%d), lon %.1f, lat %.1f",i,j,din$longitude[i],din$latitude[j]))
      #    plot.tser.cubic(x=x.time.ext,ybar=gpp.cell.ext,fit=fit.gpp,main=sprintf("GPP at ij=(%d,%d), lon %.1f, lat %.1f",i,j,din$lon[i],din$lat[j]))

      #    plot.tser.cubic(x=x.time.ext,ybar=rtot.cell.ext,fit=fit.rtot,main=sprintf("Rtot at ij=(%d,%d), lon %.1f, lat %.1f",i,j,din$lon[i],din$lat[j]))
      dev.off()
    }

    # Strip pad coefficients before storing; output dims stay c(360,180,nmon).
    piqsfit.gpp$a[i,j,] <- fit.gpp$a[keep]
    piqsfit.gpp$b[i,j,] <- fit.gpp$b[keep]
    piqsfit.gpp$c[i,j,] <- fit.gpp$c[keep]
#    piqsfit.gpp$d[i,j,] <- fit.gpp$d[keep]

    piqsfit.resp$a[i,j,] <- fit.rtot$a[keep]
    piqsfit.resp$b[i,j,] <- fit.rtot$b[keep]
    piqsfit.resp$c[i,j,] <- fit.rtot$c[keep]
#    piqsfit.resp$d[i,j,] <- fit.rtot$d[keep]
    
    pb <- progress.bar.print(pb,ipb)
    
  }
}

progress.bar.end(pb)

piqsfit.time <- x.time[1:(length(x.time)-1)]

# Metadata so downstream consumers (diurnalize-ERA5.r and any future readers)
# can tell what padding was applied. Older .rda files without this object
# should be treated as pad.left = pad.right = 0.
piqsfit.meta <- list(pad.left  = pad.left,
                     pad.right = pad.right,
                     fit.range = range(x.time.ext),
                     saved.range = range(piqsfit.time),
                     written.at = format(Sys.time(), tz="UTC", usetz=TRUE))

save(file='fit.piqs.rda',piqsfit.gpp,piqsfit.resp,piqsfit.time,piqsfit.meta)
