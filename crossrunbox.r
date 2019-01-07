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

# plot of sensitivities and specificities, LR+ and LR-:
# plotting parameters:
col.0 <- "black"
col.4 <- "brown"
col.8 <- "blue"
col1.2 <- "red"
col.targ <- "violet"
y.down <- .2
y.targ <- 1
x0 <- 0
x1 <- 20

# the plots:
par(mfrow=c(2,2))
par(mar=c(bottom=3, left=4, top=0, right=0)+.1)
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(0,1), xlab="", ylab="")
lines(x=c(x0,x1), y=c(y.targ,y.targ), col=col.targ)
text(x=x1, y=y.targ, col=col.targ, pos=4, labels="target specificity (0.925)")
lines(x=c(x0,x1), y=c(y.targ-y.down,y.targ-y.down), col=col.0)
text(x=x1, y=y.targ-y.down, col=col.0, pos=4, labels="specificity (dashed: 'best' box)")
lines(x=c(x0,x1), y=c(y.targ-2*y.down,y.targ-2*y.down), col=col1.2)
text(x=x1, y=y.targ-2*y.down, col=col1.2, pos=4, labels="sensitivity and LR,\nshift 1.2 (dashed: 'best' box)")
lines(x=c(x0,x1), y=c(y.targ-3*y.down,y.targ-3*y.down), col=col.8)
text(x=x1, y=y.targ-3*y.down, col=col.8, pos=4, 
     labels="sensitivity and LR,\n(target) shift 0.8 (dashed: 'best' box)")
lines(x=c(x0,x1), y=c(y.targ-4*y.down,y.targ-4*y.down), col=col.4)
text(x=x1, y=y.targ-4*y.down, col=col.4, pos=4, labels="sensitivity and LR,\nshift 0.4 (dashed: 'best' box)")

plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(0,1), xlab="", ylab="specificity and sensitivites")
box()
axis(1, at=10*c(1:10))
axis(2, las=1., at=c(0,.2,.4,.6,.8,.925,1),
     labels=c("0","0.2","0.4","0.6","0.8","0.925","1"))
mtext(side=1,line=2.2, text="length of series")
abline(h=.925, col=col.targ)
points(x=10:100, y=bounds$pa0, type="l", col=col.0)
points(x=10:100, y=bounds$pb0, type="l", lty="dashed", col=col.0)
points(x=10:100, y=1-bounds$pa.4, type="l", col=col.4)
points(x=10:100, y=1-bounds$pb.4, type="l", lty="dashed", col=col.4)
points(x=10:100, y=1-bounds$pa.8, type="l", col=col.8)
points(x=10:100, y=1-bounds$pb.8, type="l", lty="dashed", col=col.8)
points(x=10:100, y=1-bounds$pa1.2, type="l", col=col1.2)
points(x=10:100, y=1-bounds$pb1.2, type="l", lty="dashed", col=col1.2)

# graph of positive likelihood ratios:
par(mar=c(bottom=3, left=4, top=0, right=0)+.1)
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(0,20), xlab="",
     ylab="Positive likelihood ratio")
points(x=10:100, y=bounds$lrposa.4, type="l", col=col.4)
points(x=10:100, y=bounds$lrposb.4, type="l", lty="dashed", col=col.4)
points(x=10:100, y=bounds$lrposa.8, type="l", col=col.8)
points(x=10:100, y=bounds$lrposb.8, type="l", lty="dashed", col=col.8)
points(x=10:100, y=bounds$lrposa1.2, type="l", col=col1.2)
points(x=10:100, y=bounds$lrposb1.2, type="l", lty="dashed", col=col1.2)
box()
axis(1)
axis(2, las=1)
par(mar=c(bottom=5, left=4, top=4, right=2)+.1)

# graph of negative likelihood ratios:
par(mar=c(bottom=3, left=4, top=0, right=0)+.1)
plot(x=10:100, rep(.5,91), axes=FALSE, col="white", 
     ylim=c(0,1), xlab="",
     ylab="Negative likelihood ratio")
points(x=10:100, y=bounds$lrnega.4, type="l", col=col.4)
points(x=10:100, y=bounds$lrnegb.4, type="l", lty="dashed", col=col.4)
points(x=10:100, y=bounds$lrnega.8, type="l", col=col.8)
points(x=10:100, y=bounds$lrnegb.8, type="l", lty="dashed", col=col.8)
points(x=10:100, y=bounds$lrnega1.2, type="l", col=col1.2)
points(x=10:100, y=bounds$lrnegb1.2, type="l", lty="dashed", col=col1.2)
box()
axis(1)
axis(2, las=1)
par(mar=c(bottom=5, left=4, top=4, right=2)+.1)
par(mfrow=c(1,1))

# save functions, function bestbox and plotting paramters:
save(bounds, bestbox, ca1, cb1, col.0, col.4, col.8, col.targ, col1.2,
     la1, lb1, x0, x1, y.down, y.targ,
     file="crossrunbox1.Rdata")
