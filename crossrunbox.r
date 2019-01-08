library(Rmpfr)
library(crossrun)

# compute simultaneous distributions for n=1, ..., 100 in the symmetric case:
Sys.time()
cr100 <- crossrunsymm(100)
Sys.time() # about 3 minutes (R3.4.4, Rmpfr 0-7.0)

# shift 0.2 to 3:
Sys.time()
crs100..2 <- crossrunshift(shift=.2, printn=TRUE)$pt
crs100..4 <- crossrunshift(shift=.4, printn=TRUE)$pt
crs100..6 <- crossrunshift(shift=.6, printn=TRUE)$pt
crs100..8 <- crossrunshift(shift=.8, printn=TRUE)$pt
crs100.1 <- crossrunshift(shift=1, printn=TRUE)$pt
crs100.1.2 <- crossrunshift(shift=1.2, printn=TRUE)$pt
crs100.1.4 <- crossrunshift(shift=1.4, printn=TRUE)$pt
crs100.1.6 <- crossrunshift(shift=1.6, printn=TRUE)$pt
crs100.1.8 <- crossrunshift(shift=1.8, printn=TRUE)$pt
crs100.2 <- crossrunshift(shift=2, printn=TRUE)$pt
crs100.2.2 <- crossrunshift(shift=2.2, printn=TRUE)$pt
crs100.2.4 <- crossrunshift(shift=2.4, printn=TRUE)$pt
crs100.2.6 <- crossrunshift(shift=2.6, printn=TRUE)$pt
crs100.2.8 <- crossrunshift(shift=2.8, printn=TRUE)$pt
crs100.3 <- crossrunshift(shift=3, printn=TRUE)$pt
Sys.time()


# Table 1 in "Run charts revisited", PLOS ONE November 25, 2014:
bounds <- data.frame(
  n=10:100, 
  ca=rep(2:41,c(2,2,3,2,3,2,3,2,2,3,2,2,3, # 41
                2,2,3,2,2,2,3,2,2,2,3, # 65
                2,2,2,3,rep(2,4),3,rep(2,4),3,2,2)),
  la=c(rep(6,2), rep(7,22-11), rep(8,45-22), rep(9,90-45), 
       rep(10,100-90)))
bounds

# function for box with lowest probability for the target shift,
# among boxes with probability >=  target for shift 0.
# only boxes with positive corner probability considered.
# subsequent deletion within the border if possible is placed
# in separate function taken c and l for this box as input:
bestbox <- function(pt0=cr100$pt, pts=crs100..8, 
                     target=.925, n1=100, mult=2, prec=120) {
  nill <- mpfr(0,prec);one <- mpfr(1,prec);two <- mpfr(2,prec)
  multm <- mpfr(mult,prec);targetm <- mpfr(target,prec)
  targt <- targetm*(multm^(n1-1)) # target on "times" scale
  pt0n <- pt0[[n1]]
  ptsn <- pts[[n1]]
  bpt0 <- boxprobt(pt0n) # box probabilities for no shift
  bpttarg <- boxprobt(ptsn) # box probabilities for target shift
  boxprt <- two*(multm^(n1-1)) # initialize to impossible high value
  for (cc in 0:(n1-1)) for (ll in 1:n1) {
    if (pt0n[cc+1,ll]>nill & bpt0[cc+1,ll]>=targt & bpttarg[cc+1,ll]<boxprt) {
      c1 <- cc
      l1 <- ll
      boxprt <- bpttarg[cc+1,ll]
    }
  } # end search through (c,l) with positive no (,l)
  return(c(c1,l1)) # only box som far, reduction later
} # end function bestbox

# find alternative ("best") boxes:
bounds$cb <- NA
bounds$lb <- NA
for (nn in 10:100) {
  print(Sys.time())
  print(nn)
  bounds[bounds$n==nn,c("cb","lb")] <- bestbox(n1=nn)
  }
bounds[,c("n","ca","cb","la","lb")]

# find box (no signal) probabilities for shift 0,.4,.8,1.2
# for original and "best" boxes:
bounds$pa0 <- NA
bounds$pa.4 <- NA
bounds$pa.8 <- NA
bounds$pa1.2 <- NA
bounds$pb0 <- NA
bounds$pb.4 <- NA
bounds$pb.8 <- NA
bounds$pb1.2 <- NA
for (nn in 10:100) {
  print(Sys.time())
  print(nn)
  ca1 <- bounds$ca[bounds$n==nn]
  la1 <- bounds$la[bounds$n==nn]
  cb1 <- bounds$cb[bounds$n==nn]
  lb1 <- bounds$lb[bounds$n==nn]
  bounds$pa0[bounds$n==nn] <- 
    as.numeric(sum(cr100$pt[[nn]][(ca1+1):nn,1:la1])/sum(cr100$pt[[nn]]))
  bounds$pa.4[bounds$n==nn] <- 
    as.numeric(sum(crs100..4[[nn]][(ca1+1):nn,1:la1])/sum(crs100..4[[nn]]))
  bounds$pa.8[bounds$n==nn] <- 
    as.numeric(sum(crs100..8[[nn]][(ca1+1):nn,1:la1])/sum(crs100..8[[nn]]))
  bounds$pa1.2[bounds$n==nn] <- 
    as.numeric(sum(crs100.1.2[[nn]][(ca1+1):nn,1:la1])/sum(crs100.1.2[[nn]]))
  bounds$pb0[bounds$n==nn] <- 
    as.numeric(sum(cr100$pt[[nn]][(cb1+1):nn,1:lb1])/sum(cr100$pt[[nn]]))
  bounds$pb.4[bounds$n==nn] <- 
    as.numeric(sum(crs100..4[[nn]][(cb1+1):nn,1:lb1])/sum(crs100..4[[nn]]))
  bounds$pb.8[bounds$n==nn] <- 
    as.numeric(sum(crs100..8[[nn]][(cb1+1):nn,1:lb1])/sum(crs100..8[[nn]]))
  bounds$pb1.2[bounds$n==nn] <- 
    as.numeric(sum(crs100.1.2[[nn]][(cb1+1):nn,1:lb1])/sum(crs100.1.2[[nn]]))
}
summary(bounds)

# computation of likelihood ratios:
bounds$lrposa.4 <- (1-bounds$pa.4)/(1-bounds$pa0)
bounds$lrnega.4 <- bounds$pa.4/bounds$pa0
bounds$lrposb.4 <- (1-bounds$pb.4)/(1-bounds$pb0)
bounds$lrnegb.4 <- bounds$pb.4/bounds$pb0
bounds$lrposa.8 <- (1-bounds$pa.8)/(1-bounds$pa0)
bounds$lrnega.8 <- bounds$pa.8/bounds$pa0
bounds$lrposb.8 <- (1-bounds$pb.8)/(1-bounds$pb0)
bounds$lrnegb.8 <- bounds$pb.8/bounds$pb0
bounds$lrposa1.2 <- (1-bounds$pa1.2)/(1-bounds$pa0)
bounds$lrnega1.2 <- bounds$pa1.2/bounds$pa0
bounds$lrposb1.2 <- (1-bounds$pb1.2)/(1-bounds$pb0)
bounds$lrnegb1.2 <- bounds$pb1.2/bounds$pb0
summary(bounds[,14:25])

# function for cutting a box while keeping probability >=
# target for shift 0. No cutting if the corner cannot be 
# removed. If the corner may be removed, it is attempted 
# to remove parts of the border, starting from the corner,
# in the direction with highest point probability for the
# target shift:
cutbox <- function(pt0=cr100$pt, pts=crs100..8, target=.925, 
                   n1=100, c1=41, l1=10, mult=2, prec=120) {
  nill <- mpfr(0,prec);one <- mpfr(1,prec);two <- mpfr(2,prec)
  multm <- mpfr(mult,prec);targetm <- mpfr(target,prec)
  targt <- targetm*(multm^(n1-1)) # target on "times" scale
  pt0n <- pt0[[n1]];ptsn <- pts[[n1]]
  bpt0 <- boxprobt(pt0n) # box probabilities for no shift, pt scale
  boxpt0 <- bpt0[c1+1,l1] # no shift probability of actual box, pt scale
  cornerpt0 <- pt0n[c1+1,l1] # no shift corner probability, pt scale
  finished <- FALSE
  cbord <- NA;lbord <- NA
  if (boxpt0-cornerpt0>=targt) {
    cutboxpt0 <- boxpt0-cornerpt0 # pt of cutted box after removed corner
    cbord <- c1+1;lbord <- l1-1
    while (!finished) {
      pt0n.directionc <- pt0n[cbord+1,l1]
      pt0n.directionl <- pt0n[c1+1,lbord] 
      ptsn.directionc <- ptsn[cbord+1,l1]
      ptsn.directionl <- ptsn[c1+1,lbord] 
      if ((cutboxpt0-pt0n.directionc<targt|pt0n.directionc==0)&
          (cutboxpt0-pt0n.directionl<targt|pt0n.directionl==0)) {
        finished <- TRUE 
      } else 
        if (cutboxpt0-pt0n.directionc<targt|pt0n.directionc==0) {
          lstrip <- pt0n[c1+1,lbord:1]
          nlstrip <- length(lstrip)
          maxlstrip <- max((1:nlstrip)[lstrip>0])
          lstrip <- lstrip[(1:nlstrip)<=maxlstrip]
          lstripcum <- cumsum(lstrip)
          if (cutboxpt0-max(lstripcum)>=targt) lbord <- 0 else # 0 cannot occurr
            lbord <- lbord+1-min((1:nlstrip)[cutboxpt0-lstripcum<targt])
          finished <- TRUE
        } else if (cutboxpt0-pt0n.directionl<targt|pt0n.directionl==0) {
          cstrip <- pt0n[(cbord+1):n1,l1]
          ncstrip <- length(cstrip)
          maxcstrip <- max((1:ncstrip)[cstrip>0])
          cstrip <- cstrip[(1:ncstrip)<=maxcstrip]
          cstripcum <- cumsum(cstrip)
          if (cutboxpt0-max(cstripcum)>=targt) cbord <- n1 else # n1 cannot occurr
            cbord <- cbord+min((1:ncstrip)[cutboxpt0-cstripcum<targt])-1
          finished <- TRUE
        } else if (ptsn.directionc>=ptsn.directionl) {
          cbord <- cbord+1
          cutboxpt0 <- cutboxpt0 - pt0n.directionc
        } else if (ptsn.directionc<ptsn.directionl) {
          lbord <- lbord-1
          cutboxpt0 <- cutboxpt0 - pt0n.directionl
        }
    } # end while loop
  } # end if corner may be removed
  return(c(cbord,lbord))
} # end function cutbox

# find cutted  boxes:
bounds$cbord <- NA
bounds$lbord <- NA
for (nn in 10:100) {
  print(Sys.time())
  print(nn)
  bounds[bounds$n==nn,c("cbord","lbord")] <- 
    cutbox(n1=nn,c1=bounds$cb[bounds$n==nn],l1=bounds$lb[bounds$n==nn])
}
data.frame(bounds[,c("ca","la","cb","lb","cbord","lbord")], row.names=bounds$n)

# find cutbox probabilities for shift 0,.4,.8,1.2:
bounds$pc0 <- NA
bounds$pc.4 <- NA
bounds$pc.8 <- NA
bounds$pc1.2 <- NA
for (nn in 10:100) {
  print(Sys.time())
  print(nn)
  cb1 <- bounds$cb[bounds$n==nn]
  lb1 <- bounds$lb[bounds$n==nn]
  cbord1 <- bounds$cbord[bounds$n==nn]
  lbord1 <- bounds$lbord[bounds$n==nn]
  if (is.na(cbord1)==1) {
    bounds$pc0[bounds$n==nn] <- bounds$pb0[bounds$n==nn]
    bounds$pc.4[bounds$n==nn] <- bounds$pb.4[bounds$n==nn]
    bounds$pc.8[bounds$n==nn] <- bounds$pb.8[bounds$n==nn]
    bounds$pc1.2[bounds$n==nn] <- bounds$pb1.2[bounds$n==nn]
  } else {
    bounds$pc0[bounds$n==nn] <- 
      as.numeric((sum(cr100$pt[[nn]][(cb1+2):nn,1:(lb1-1)]) +
                    sum(cr100$pt[[nn]][(cbord1+1):nn,lb1]) +
                    sum(cr100$pt[[nn]][cb1+1,1:lbord1]))/
                   sum(cr100$pt[[nn]]))
    bounds$pc.4[bounds$n==nn] <- 
      as.numeric((sum(crs100..4[[nn]][(cb1+2):nn,1:(lb1-1)]) +
                    sum(crs100..4[[nn]][(cbord1+1):nn,lb1]) +
                    sum(crs100..4[[nn]][cb1+1,1:lbord1]))/
                   sum(crs100..4[[nn]]))
    bounds$pc.8[bounds$n==nn] <- 
      as.numeric((sum(crs100..8[[nn]][(cb1+2):nn,1:(lb1-1)]) +
                    sum(crs100..8[[nn]][(cbord1+1):nn,lb1]) +
                    sum(crs100..8[[nn]][cb1+1,1:lbord1]))/
                   sum(crs100..8[[nn]]))
    bounds$pc1.2[bounds$n==nn] <- 
      as.numeric((sum(crs100.1.2[[nn]][(cb1+2):nn,1:(lb1-1)]) +
                    sum(crs100.1.2[[nn]][(cbord1+1):nn,lb1]) +
                    sum(crs100.1.2[[nn]][cb1+1,1:lbord1]))/
                   sum(crs100.1.2[[nn]]))
  }}
summary(bounds)

# computation of likelihood ratios for cutboxes:
bounds$lrposc.4 <- (1-bounds$pc.4)/(1-bounds$pc0)
bounds$lrnegc.4 <- bounds$pc.4/bounds$pc0
bounds$lrposc.8 <- (1-bounds$pc.8)/(1-bounds$pc0)
bounds$lrnegc.8 <- bounds$pc.8/bounds$pc0
bounds$lrposc1.2 <- (1-bounds$pc1.2)/(1-bounds$pc0)
bounds$lrnegc1.2 <- bounds$pc1.2/bounds$pc0
summary(bounds[,c("lrposa.4","lrposb.4","lrposc.4","lrposa.8","lrposb.8","lrposc.8",
                  "lrposa1.2","lrposb1.2","lrposc1.2",
                  "lrnega.4","lrnegb.4","lrnegc.4","lrnega.8","lrnegb.8","lrnegc.8",
                  "lrnega1.2","lrnegb1.2","lrnegc1.2")])
# sequence lengths where original rules are below target:
below.target <- c(10:100)[bounds$pa0<.925]
n.below.target <- length(below.target)

# plots:
cola <- "red"
colb <- "orange"
colc <- "blue"
# specificity:
par(mar=c(bottom=3, left=4, top=0, right=0)+.1)
uppspec <- ceiling(max(c(bounds$pa0,bounds$pb0,
                             bounds$pc0))*100)/100
lowspec <- floor(min(c(bounds$pa0,bounds$pb0,
                           bounds$pc0))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowspec,uppspec), xlab="", ylab="specificity")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
abline(h=.925, col=col.targ)
points(x=10:100, y=bounds$pa0, type="l", col=cola)
points(x=10:100, y=bounds$pb0, type="l", col=colb)
points(x=10:100, y=bounds$pc0, type="l", col=colc)
points(x=below.target, y=rep(.925,n.below.target), pch=19)
# sensitivity p=0.4:
par(mar=c(bottom=2, left=4, top=0, right=0)+.1)
uppsens.4 <- ceiling(max(1-c(bounds$pa.4,bounds$pb.4,
                             bounds$pc.4))*100)/100
lowsens.4 <- floor(min(1-c(bounds$pa.4,bounds$pb.4,
                           bounds$pc.4))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowsens.4,uppsens.4), xlab="", 
     ylab="sensitivity shift 0.4")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
points(x=10:100, y=1-bounds$pa.4, type="l", col=cola)
points(x=10:100, y=1-bounds$pb.4, type="l", col=colb)
points(x=10:100, y=1-bounds$pc.4, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(1-bounds$pa.4,1-bounds$pb.4,1-bounds$pc.4)[below.target-9])
# sensitivity p=0.8:
uppsens.8 <- ceiling(max(1-c(bounds$pa.8,bounds$pb.8,
                             bounds$pc.8))*100)/100
lowsens.8 <- floor(min(1-c(bounds$pa.8,bounds$pb.8,
                           bounds$pc.8))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowsens.8,uppsens.8), xlab="", 
     ylab="sensitivity shift 0.8")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=1-bounds$pa.8, type="l", col=cola)
points(x=10:100, y=1-bounds$pb.8, type="l", col=colb)
points(x=10:100, y=1-bounds$pc.8, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(1-bounds$pa.8,1-bounds$pb.8,1-bounds$pc.8)[below.target-9])
# sensitivity p=1.2:
uppsens1.2 <- ceiling(max(1-c(bounds$pa1.2,bounds$pb1.2,
                             bounds$pc1.2))*100)/100
lowsens1.2 <- floor(min(1-c(bounds$pa1.2,bounds$pb1.2,
                           bounds$pc1.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowsens1.2,uppsens1.2), xlab="", 
     ylab="sensitivity shift 1.2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=1-bounds$pa1.2, type="l", col=cola)
points(x=10:100, y=1-bounds$pb1.2, type="l", col=colb)
points(x=10:100, y=1-bounds$pc1.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(1-bounds$pa1.2,1-bounds$pb1.2,1-bounds$pc1.2)[below.target-9])
# positive likelihood ratios p=0.4:
upplrpos.4 <- ceiling(max(c(bounds$lrposa.4,bounds$lrposb.4,
                            bounds$lrposc.4))*100)/100
lowlrpos.4 <- floor(min(c(bounds$lrposa.4,bounds$lrposb.4,
                          bounds$lrposc.4))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrpos.4,upplrpos.4), xlab="", 
     ylab="positive lhood ratio shift 0.4")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrposa.4, type="l", col=cola)
points(x=10:100, y=bounds$lrposb.4, type="l", col=colb)
points(x=10:100, y=bounds$lrposc.4, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrposa.4,bounds$lrposb.4,bounds$lrposc.4)[below.target-9])
# positive likelihood ratios p=0.8:
upplrpos.8 <- ceiling(max(c(bounds$lrposa.8,bounds$lrposb.8,
                            bounds$lrposc.8))*100)/100
lowlrpos.8 <- floor(min(c(bounds$lrposa.8,bounds$lrposb.8,
                          bounds$lrposc.8))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrpos.8,upplrpos.8), xlab="", 
     ylab="positive lhood ratio shift 0.8")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrposa.8, type="l", col=cola)
points(x=10:100, y=bounds$lrposb.8, type="l", col=colb)
points(x=10:100, y=bounds$lrposc.8, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrposa.8,bounds$lrposb.8,bounds$lrposc.8)[below.target-9])
# positive likelihood ratios p=1.2:
upplrpos1.2 <- ceiling(max(c(bounds$lrposa1.2,bounds$lrposb1.2,
                            bounds$lrposc1.2))*100)/100
lowlrpos1.2 <- floor(min(c(bounds$lrposa1.2,bounds$lrposb1.2,
                          bounds$lrposc1.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrpos1.2,upplrpos1.2), xlab="", 
     ylab="positive lhood ratio shift 1.2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrposa1.2, type="l", col=cola)
points(x=10:100, y=bounds$lrposb1.2, type="l", col=colb)
points(x=10:100, y=bounds$lrposc1.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrposa1.2,bounds$lrposb1.2,bounds$lrposc1.2)[below.target-9])
# negative likelihood ratios p=0.4:
upplrneg.4 <- ceiling(max(c(bounds$lrnega.4,bounds$lrnegb.4,
                            bounds$lrnegc.4))*100)/100
lowlrneg.4 <- floor(min(c(bounds$lrnega.4,bounds$lrnegb.4,
                          bounds$lrnegc.4))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrneg.4,upplrneg.4), xlab="", 
     ylab="negative lhood ratio shift 0.4")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrnega.4, type="l", col=cola)
points(x=10:100, y=bounds$lrnegb.4, type="l", col=colb)
points(x=10:100, y=bounds$lrnegc.4, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrnega.4,bounds$lrnegb.4,bounds$lrnegc.4)[below.target-9])
# negative likelihood ratios p=0.8:
upplrneg.8 <- ceiling(max(c(bounds$lrnega.8,bounds$lrnegb.8,
                            bounds$lrnegc.8))*100)/100
lowlrneg.8 <- floor(min(c(bounds$lrnega.8,bounds$lrnegb.8,
                          bounds$lrnegc.8))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrneg.8,upplrneg.8), xlab="", 
     ylab="negative lhood ratio shift 0.8")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrnega.8, type="l", col=cola)
points(x=10:100, y=bounds$lrnegb.8, type="l", col=colb)
points(x=10:100, y=bounds$lrnegc.8, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrnega.8,bounds$lrnegb.8,bounds$lrnegc.8)[below.target-9])
# negative likelihood ratios p=1.2:
upplrneg1.2 <- ceiling(max(c(bounds$lrnega1.2,bounds$lrnegb1.2,
                             bounds$lrnegc1.2))*100)/100
lowlrneg1.2 <- floor(min(c(bounds$lrnega1.2,bounds$lrnegb1.2,
                           bounds$lrnegc1.2))*100)/100
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(lowlrneg1.2,upplrneg1.2), xlab="", 
     ylab="negative lhood ratio shift 1.2")
box()
axis(1, at=10*c(1:10))
axis(2, las=1)
mtext(side=1,line=2.2, text="length of sequence")
points(x=10:100, y=bounds$lrnega1.2, type="l", col=cola)
points(x=10:100, y=bounds$lrnegb1.2, type="l", col=colb)
points(x=10:100, y=bounds$lrnegc1.2, type="l", col=colc)
points(x=below.target, pch=19,
       y=pmax(bounds$lrnega1.2,bounds$lrnegb1.2,bounds$lrnegc1.2)[below.target-9])

# save all except large objects:
save(list=ls()[substr(ls(),1,2)!="cr"], file="crossrunbox1.Rdata")

