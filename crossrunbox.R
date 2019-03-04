library(Rmpfr)
library(crossrun)

# function for box with lowest probability for the target shift, among boxes
# with probability >=  target for shift 0. only boxes with positive corner
# probability considered. subsequent deletion within the border if possible is
# placed in separate function taken c and l for this box as input:
bestbox <- function(pt0    = cr100,
                    pts    = crs100.0.8,
                    target = 0.925,
                    n1     = 100,
                    mult   = 2,
                    prec   = 120) {
  nill    <- mpfr(0, prec)
  one     <- mpfr(1, prec)
  two     <- mpfr(2, prec)
  multm   <- mpfr(mult, prec)
  targetm <- mpfr(target, prec)
  targt   <- targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n    <- pt0[[n1]]
  ptsn    <- pts[[n1]]
  bpt0    <- boxprobt(pt0n) # box probabilities for no shift
  bpttarg <- boxprobt(ptsn) # box probabilities for target shift
  boxprt  <- two * (multm ^ (n1 - 1)) # initialize to impossible high value
  
  for (cc in 0:(n1 - 1))
    for (ll in 1:n1) {
      if (pt0n[cc + 1, ll] > nill &
          bpt0[cc + 1, ll] >= targt & bpttarg[cc + 1, ll] < boxprt) {
        c1 <- cc
        l1 <- ll
        boxprt <- bpttarg[cc + 1, ll]
      }
    } # end search through (c,l) with positive no (,l)
  return(c(c1, l1)) # only box som far, reduction later
} # end function bestbox

# function for cutting a box while keeping probability >= target for shift 0. No
# cutting if the corner cannot be removed. If the corner may be removed, it is
# attempted to remove parts of the border, starting from the corner, in the
# direction with highest point probability for the target shift:
cutbox <- function(pt0    = cr100,
                   pts    = crs100.0.8,
                   target = 0.925,
                   n1     = 100,
                   c1     = 41,
                   l1     = 10,
                   mult   = 2,
                   prec   = 120) {
  nill      <- mpfr(0, prec)
  one       <- mpfr(1, prec)
  two       <- mpfr(2, prec)
  multm     <- mpfr(mult, prec)
  targetm   <- mpfr(target, prec)
  targt     <- targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n      <- pt0[[n1]]
  ptsn      <- pts[[n1]]
  bpt0      <- boxprobt(pt0n) # box probabilities for no shift, pt scale
  boxpt0    <- bpt0[c1 + 1, l1] # no shift probability of actual box, pt scale
  cornerpt0 <- pt0n[c1 + 1, l1] # no shift corner probability, pt scale
  finished  <- FALSE
  cbord     <- NA
  lbord     <- NA
  
  if (boxpt0 - cornerpt0 >= targt) {
    cutboxpt0 <- boxpt0 - cornerpt0 # pt of cutted box after removed corner
    cbord     <- c1 + 1
    lbord     <- l1 - 1
    
    while (!finished) {
      pt0n.directionc <- pt0n[cbord + 1, l1]
      pt0n.directionl <- pt0n[c1 + 1, lbord]
      ptsn.directionc <- ptsn[cbord + 1, l1]
      ptsn.directionl <- ptsn[c1 + 1, lbord]
      if ((cutboxpt0 - pt0n.directionc < targt | pt0n.directionc == 0) &
          (cutboxpt0 - pt0n.directionl < targt |
           pt0n.directionl == 0)) {
        finished <- TRUE
      } else if (cutboxpt0 - pt0n.directionc < targt | pt0n.directionc == 0) {
        lstrip    <- pt0n[c1 + 1, lbord:1]
        nlstrip   <- length(lstrip)
        maxlstrip <- max((1:nlstrip)[lstrip > 0])
        lstrip    <- lstrip[(1:nlstrip) <= maxlstrip]
        lstripcum <- cumsum(lstrip)
        if (cutboxpt0 - max(lstripcum) >= targt){
          lbord <- 0
        } else {
          # 0 cannot occurr
          lbord <- lbord + 1 - min((1:nlstrip)[cutboxpt0 - lstripcum < targt])
        }
        finished <- TRUE
      } else if (cutboxpt0 - pt0n.directionl < targt | pt0n.directionl == 0) {
        cstrip    <- pt0n[(cbord + 1):n1, l1]
        ncstrip   <- length(cstrip)
        maxcstrip <- max((1:ncstrip)[cstrip > 0])
        cstrip    <- cstrip[(1:ncstrip) <= maxcstrip]
        cstripcum <- cumsum(cstrip)
        if (cutboxpt0 - max(cstripcum) >= targt) {
          cbord <- n1
        } else {
          # n1 cannot occurr
          cbord <- cbord + min((1:ncstrip)[cutboxpt0 - cstripcum < targt]) - 1
        }
        finished <- TRUE
      } else if (ptsn.directionc >= ptsn.directionl) {
        cbord     <- cbord + 1
        cutboxpt0 <- cutboxpt0 - pt0n.directionc
      } else if (ptsn.directionc < ptsn.directionl) {
        lbord     <- lbord - 1
        cutboxpt0 <- cutboxpt0 - pt0n.directionl
      }
    } # end while loop
  } # end if corner may be removed
  return(c(cbord, lbord))
} # end function cutbox

# compute simultaneous distributions for n=1, ..., 100 in the symmetric case:
cr100 <- crossrunsymm(100, printn = TRUE)$pt

# shift 0.2 to 3:
crs100.0.2 <- crossrunshift(shift = 0.2, printn = TRUE)$pt
crs100.0.4 <- crossrunshift(shift = 0.4, printn = TRUE)$pt
crs100.0.6 <- crossrunshift(shift = 0.6, printn = TRUE)$pt
crs100.0.8 <- crossrunshift(shift = 0.8, printn = TRUE)$pt
crs100.1.0 <- crossrunshift(shift = 1.0, printn = TRUE)$pt
crs100.1.2 <- crossrunshift(shift = 1.2, printn = TRUE)$pt
crs100.1.4 <- crossrunshift(shift = 1.4, printn = TRUE)$pt
crs100.1.6 <- crossrunshift(shift = 1.6, printn = TRUE)$pt
crs100.1.8 <- crossrunshift(shift = 1.8, printn = TRUE)$pt
crs100.2.0 <- crossrunshift(shift = 2.0, printn = TRUE)$pt
crs100.2.2 <- crossrunshift(shift = 2.2, printn = TRUE)$pt
crs100.2.4 <- crossrunshift(shift = 2.4, printn = TRUE)$pt
crs100.2.6 <- crossrunshift(shift = 2.6, printn = TRUE)$pt
crs100.2.8 <- crossrunshift(shift = 2.8, printn = TRUE)$pt
crs100.3.0 <- crossrunshift(shift = 3.0, printn = TRUE)$pt

# Table 1 in "Run charts revisited", PLOS ONE November 25, 2014:
bounds <- data.frame(
  n = 10:100,
  ca = qbinom(0.05, 10:100 - 1, 0.5),
  la = round(log2(10:100) + 3))
row.names(bounds) <- bounds$n
# bounds

# find alternative ("best") boxes:
bounds$cb <- NA
bounds$lb <- NA

for (nn in 10:100) {
  print(nn)
  
  bounds[bounds$n == nn, c("cb", "lb")] <- bestbox(n1 = nn)
}

# bounds[, c("n","ca","cb","la","lb")]

# find cutted  boxes:
bounds$cbord <- NA
bounds$lbord <- NA

for (nn in 10:100) {
  print(nn)
  
  bounds[bounds$n == nn, c("cbord", "lbord")] <-
    cutbox(n1 = nn,
           c1 = bounds$cb[bounds$n == nn],
           l1 = bounds$lb[bounds$n == nn])
}

# bounds

# find no signal probabilities for the shifts, original and besbox rules:
bounds$pa_0.0 <- NA
bounds$pa_0.2 <- NA
bounds$pa_0.4 <- NA
bounds$pa_0.6 <- NA
bounds$pa_0.8 <- NA
bounds$pa_1.0 <- NA
bounds$pa_1.2 <- NA
bounds$pa_1.4 <- NA
bounds$pa_1.6 <- NA
bounds$pa_1.8 <- NA
bounds$pa_2.0 <- NA
bounds$pa_2.2 <- NA
bounds$pa_2.4 <- NA
bounds$pa_2.6 <- NA
bounds$pa_2.8 <- NA
bounds$pa_3.0 <- NA
bounds$pb_0.0 <- NA
bounds$pb_0.2 <- NA
bounds$pb_0.4 <- NA
bounds$pb_0.6 <- NA
bounds$pb_0.8 <- NA
bounds$pb_1.0 <- NA
bounds$pb_1.2 <- NA
bounds$pb_1.4 <- NA
bounds$pb_1.6 <- NA
bounds$pb_1.8 <- NA
bounds$pb_2.0 <- NA
bounds$pb_2.2 <- NA
bounds$pb_2.4 <- NA
bounds$pb_2.6 <- NA
bounds$pb_2.8 <- NA
bounds$pb_3.0 <- NA

for (nn in 10:100) {
  print(nn)
  
  ca1 <- bounds$ca[bounds$n == nn]
  la1 <- bounds$la[bounds$n == nn]
  cb1 <- bounds$cb[bounds$n == nn]
  lb1 <- bounds$lb[bounds$n == nn]
  bounds$pa_0.0[bounds$n == nn] <-
    as.numeric(sum(cr100$pt[[nn]][(ca1 + 1):nn, 1:la1]) / sum(cr100$pt[[nn]]))
  bounds$pa_0.2[bounds$n == nn] <-
    as.numeric(sum(crs100.0.2[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.0.2[[nn]]))
  bounds$pa_0.4[bounds$n == nn] <-
    as.numeric(sum(crs100.0.4[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.0.4[[nn]]))
  bounds$pa_0.6[bounds$n == nn] <-
    as.numeric(sum(crs100.0.6[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.0.6[[nn]]))
  bounds$pa_0.8[bounds$n == nn] <-
    as.numeric(sum(crs100.0.8[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.0.8[[nn]]))
  bounds$pa_1.0[bounds$n == nn] <-
    as.numeric(sum(crs100.1.0[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.1.0[[nn]]))
  bounds$pa_1.2[bounds$n == nn] <-
    as.numeric(sum(crs100.1.2[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.1.2[[nn]]))
  bounds$pa_1.4[bounds$n == nn] <-
    as.numeric(sum(crs100.1.4[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.1.4[[nn]]))
  bounds$pa_1.6[bounds$n == nn] <-
    as.numeric(sum(crs100.1.6[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.1.6[[nn]]))
  bounds$pa_1.8[bounds$n == nn] <-
    as.numeric(sum(crs100.1.8[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.1.8[[nn]]))
  bounds$pa_2.0[bounds$n == nn] <-
    as.numeric(sum(crs100.2.0[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.2.0[[nn]]))
  bounds$pa_2.2[bounds$n == nn] <-
    as.numeric(sum(crs100.2.2[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.2.2[[nn]]))
  bounds$pa_2.4[bounds$n == nn] <-
    as.numeric(sum(crs100.2.4[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.2.4[[nn]]))
  bounds$pa_2.6[bounds$n == nn] <-
    as.numeric(sum(crs100.2.6[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.2.6[[nn]]))
  bounds$pa_2.8[bounds$n == nn] <-
    as.numeric(sum(crs100.2.8[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.2.8[[nn]]))
  bounds$pa_3.0[bounds$n == nn] <-
    as.numeric(sum(crs100.3.0[[nn]][(ca1 + 1):nn, 1:la1]) / sum(crs100.3.0[[nn]]))
  bounds$pb_0.0[bounds$n == nn] <-
    as.numeric(sum(cr100$pt[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(cr100$pt[[nn]]))
  bounds$pb_0.2[bounds$n == nn] <-
    as.numeric(sum(crs100.0.2[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.0.2[[nn]]))
  bounds$pb_0.4[bounds$n == nn] <-
    as.numeric(sum(crs100.0.4[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.0.4[[nn]]))
  bounds$pb_0.6[bounds$n == nn] <-
    as.numeric(sum(crs100.0.6[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.0.6[[nn]]))
  bounds$pb_0.8[bounds$n == nn] <-
    as.numeric(sum(crs100.0.8[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.0.8[[nn]]))
  bounds$pb_1.0[bounds$n == nn] <-
    as.numeric(sum(crs100.1.0[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.1.0[[nn]]))
  bounds$pb_1.2[bounds$n == nn] <-
    as.numeric(sum(crs100.1.2[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.1.2[[nn]]))
  bounds$pb_1.4[bounds$n == nn] <-
    as.numeric(sum(crs100.1.4[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.1.4[[nn]]))
  bounds$pb_1.6[bounds$n == nn] <-
    as.numeric(sum(crs100.1.6[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.1.6[[nn]]))
  bounds$pb_1.8[bounds$n == nn] <-
    as.numeric(sum(crs100.1.8[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.1.8[[nn]]))
  bounds$pb_2.0[bounds$n == nn] <-
    as.numeric(sum(crs100.2.0[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.2.0[[nn]]))
  bounds$pb_2.2[bounds$n == nn] <-
    as.numeric(sum(crs100.2.2[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.2.2[[nn]]))
  bounds$pb_2.4[bounds$n == nn] <-
    as.numeric(sum(crs100.2.4[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.2.4[[nn]]))
  bounds$pb_2.6[bounds$n == nn] <-
    as.numeric(sum(crs100.2.6[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.2.6[[nn]]))
  bounds$pb_2.8[bounds$n == nn] <-
    as.numeric(sum(crs100.2.8[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.2.8[[nn]]))
  bounds$pb_3.0[bounds$n == nn] <-
    as.numeric(sum(crs100.3.0[[nn]][(cb1 + 1):nn, 1:lb1]) / sum(crs100.3.0[[nn]]))
}
# summary(bounds)

# find cutbox probabilities for the shifts:
bounds$pc_0.0 <- NA
bounds$pc_0.2 <- NA
bounds$pc_0.4 <- NA
bounds$pc_0.6 <- NA
bounds$pc_0.8 <- NA
bounds$pc_1.0 <- NA
bounds$pc_1.2 <- NA
bounds$pc_1.4 <- NA
bounds$pc_1.6 <- NA
bounds$pc_1.8 <- NA
bounds$pc_2.0 <- NA
bounds$pc_2.2 <- NA
bounds$pc_2.4 <- NA
bounds$pc_2.6 <- NA
bounds$pc_2.8 <- NA
bounds$pc_3.0 <- NA

for (nn in 10:100) {
  print(nn)
  
  cb1    <- bounds$cb[bounds$n == nn]
  lb1    <- bounds$lb[bounds$n == nn]
  cbord1 <- bounds$cbord[bounds$n == nn]
  lbord1 <- bounds$lbord[bounds$n == nn]
  
  if (is.na(cbord1) == 1) {
    bounds$pc_0.0[bounds$n == nn] <- bounds$pb_0.0[bounds$n == nn]
    bounds$pc_0.2[bounds$n == nn] <- bounds$pb_0.2[bounds$n == nn]
    bounds$pc_0.4[bounds$n == nn] <- bounds$pb_0.4[bounds$n == nn]
    bounds$pc_0.6[bounds$n == nn] <- bounds$pb_0.6[bounds$n == nn]
    bounds$pc_0.8[bounds$n == nn] <- bounds$pb_0.8[bounds$n == nn]
    bounds$pc_1.0[bounds$n == nn] <- bounds$pb_1.0[bounds$n == nn]
    bounds$pc_1.2[bounds$n == nn] <- bounds$pb_1.2[bounds$n == nn]
    bounds$pc_1.4[bounds$n == nn] <- bounds$pb_1.4[bounds$n == nn]
    bounds$pc_1.6[bounds$n == nn] <- bounds$pb_1.6[bounds$n == nn]
    bounds$pc_1.8[bounds$n == nn] <- bounds$pb_1.8[bounds$n == nn]
    bounds$pc_2.0[bounds$n == nn] <- bounds$pb_2.0[bounds$n == nn]
    bounds$pc_2.2[bounds$n == nn] <- bounds$pb_2.2[bounds$n == nn]
    bounds$pc_2.4[bounds$n == nn] <- bounds$pb_2.4[bounds$n == nn]
    bounds$pc_2.6[bounds$n == nn] <- bounds$pb_2.6[bounds$n == nn]
    bounds$pc_2.8[bounds$n == nn] <- bounds$pb_2.8[bounds$n == nn]
    bounds$pc_3.0[bounds$n == nn] <- bounds$pb_3.0[bounds$n == nn]
  } else {
    bounds$pc_0.0[bounds$n == nn] <- 
      as.numeric((sum(cr100$pt[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(cr100$pt[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(cr100$pt[[nn]][cb1 + 1, 1:lbord1]))/
                   sum(cr100$pt[[nn]]))
    bounds$pc_0.2[bounds$n == nn] <- 
      as.numeric((sum(crs100.0.2[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.0.2[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.0.2[[nn]][cb1 + 1, 1:lbord1]))/
                   sum(crs100.0.2[[nn]]))
    bounds$pc_0.4[bounds$n == nn] <- 
      as.numeric((sum(crs100.0.4[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.0.4[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.0.4[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.0.4[[nn]]))
    bounds$pc_0.6[bounds$n == nn] <-
      as.numeric((sum(crs100.0.6[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.0.6[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.0.6[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.0.6[[nn]]))
    bounds$pc_0.8[bounds$n==nn] <- 
      as.numeric((sum(crs100.0.8[[nn]][(cb1+2):nn,1:(lb1-1)]) +
                    sum(crs100.0.8[[nn]][(cbord1+1):nn,lb1]) +
                    sum(crs100.0.8[[nn]][cb1+1,1:lbord1]))/
                   sum(crs100.0.8[[nn]]))
    bounds$pc_1.0[bounds$n == nn] <-
      as.numeric((sum(crs100.1.0[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.1.0[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.1.0[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.1.0[[nn]]))
    bounds$pc_1.2[bounds$n == nn] <-
      as.numeric((sum(crs100.1.2[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.1.2[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.1.2[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.1.2[[nn]]))
    bounds$pc_1.4[bounds$n == nn] <-
      as.numeric((sum(crs100.1.4[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.1.4[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.1.4[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.1.4[[nn]]))
    bounds$pc_1.6[bounds$n == nn] <-
      as.numeric((sum(crs100.1.6[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.1.6[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.1.6[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.1.6[[nn]]))
    bounds$pc_1.8[bounds$n == nn] <-
      as.numeric((sum(crs100.1.8[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.1.8[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.1.8[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.1.8[[nn]]))
    bounds$pc_2.0[bounds$n == nn] <-
      as.numeric((sum(crs100.2.0[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.2.0[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.2.0[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.2.0[[nn]]))
    bounds$pc_2.2[bounds$n == nn] <-
      as.numeric((sum(crs100.2.2[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.2.2[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.2.2[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.2.2[[nn]]))
    bounds$pc_2.4[bounds$n == nn] <-
      as.numeric((sum(crs100.2.4[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.2.4[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.2.4[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.2.4[[nn]]))
    bounds$pc_2.6[bounds$n == nn] <-
      as.numeric((sum(crs100.2.6[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.2.6[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.2.6[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.2.6[[nn]]))
    bounds$pc_2.8[bounds$n == nn] <-
      as.numeric((sum(crs100.2.8[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.2.8[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.2.8[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.2.8[[nn]]))
    bounds$pc_3.0[bounds$n == nn] <-
      as.numeric((sum(crs100.3.0[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
                    sum(crs100.3.0[[nn]][(cbord1 + 1):nn, lb1]) +
                    sum(crs100.3.0[[nn]][cb1 + 1, 1:lbord1])) /
                   sum(crs100.3.0[[nn]]))
  }}

# computation of likelihood ratios:
bounds$lrposa_0.2 <- (1 - bounds$pa_0.2) / (1 - bounds$pa_0.0)
bounds$lrposa_0.4 <- (1 - bounds$pa_0.4) / (1 - bounds$pa_0.0)
bounds$lrposa_0.6 <- (1 - bounds$pa_0.6) / (1 - bounds$pa_0.0)
bounds$lrposa_0.8 <- (1 - bounds$pa_0.8) / (1 - bounds$pa_0.0)
bounds$lrposa_1.0 <- (1 - bounds$pa_1.0) / (1 - bounds$pa_0.0)
bounds$lrposa_1.2 <- (1 - bounds$pa_1.2) / (1 - bounds$pa_0.0)
bounds$lrposa_1.4 <- (1 - bounds$pa_1.4) / (1 - bounds$pa_0.0)
bounds$lrposa_1.6 <- (1 - bounds$pa_1.6) / (1 - bounds$pa_0.0)
bounds$lrposa_1.8 <- (1 - bounds$pa_1.8) / (1 - bounds$pa_0.0)
bounds$lrposa_2.0 <- (1 - bounds$pa_2.0) / (1 - bounds$pa_0.0)
bounds$lrposa_2.2 <- (1 - bounds$pa_2.2) / (1 - bounds$pa_0.0)
bounds$lrposa_2.4 <- (1 - bounds$pa_2.4) / (1 - bounds$pa_0.0)
bounds$lrposa_2.6 <- (1 - bounds$pa_2.6) / (1 - bounds$pa_0.0)
bounds$lrposa_2.8 <- (1 - bounds$pa_2.8) / (1 - bounds$pa_0.0)
bounds$lrposa_3.0 <- (1 - bounds$pa_3.0) / (1 - bounds$pa_0.0)
bounds$lrposb_0.2 <- (1 - bounds$pb_0.2) / (1 - bounds$pb_0.0)
bounds$lrposb_0.4 <- (1 - bounds$pb_0.4) / (1 - bounds$pb_0.0)
bounds$lrposb_0.6 <- (1 - bounds$pb_0.6) / (1 - bounds$pb_0.0)
bounds$lrposb_0.8 <- (1 - bounds$pb_0.8) / (1 - bounds$pb_0.0)
bounds$lrposb_1.0 <- (1 - bounds$pb_1.0) / (1 - bounds$pb_0.0)
bounds$lrposb_1.2 <- (1 - bounds$pb_1.2) / (1 - bounds$pb_0.0)
bounds$lrposb_1.4 <- (1 - bounds$pb_1.4) / (1 - bounds$pb_0.0)
bounds$lrposb_1.6 <- (1 - bounds$pb_1.6) / (1 - bounds$pb_0.0)
bounds$lrposb_1.8 <- (1 - bounds$pb_1.8) / (1 - bounds$pb_0.0)
bounds$lrposb_2.0 <- (1 - bounds$pb_2.0) / (1 - bounds$pb_0.0)
bounds$lrposb_2.2 <- (1 - bounds$pb_2.2) / (1 - bounds$pb_0.0)
bounds$lrposb_2.4 <- (1 - bounds$pb_2.4) / (1 - bounds$pb_0.0)
bounds$lrposb_2.6 <- (1 - bounds$pb_2.6) / (1 - bounds$pb_0.0)
bounds$lrposb_2.8 <- (1 - bounds$pb_2.8) / (1 - bounds$pb_0.0)
bounds$lrposb_3.0 <- (1 - bounds$pb_3.0) / (1 - bounds$pb_0.0)
bounds$lrposc_0.2 <- (1 - bounds$pc_0.2) / (1 - bounds$pc_0.0)
bounds$lrposc_0.4 <- (1 - bounds$pc_0.4) / (1 - bounds$pc_0.0)
bounds$lrposc_0.6 <- (1 - bounds$pc_0.6) / (1 - bounds$pc_0.0)
bounds$lrposc_0.8 <- (1 - bounds$pc_0.8) / (1 - bounds$pc_0.0)
bounds$lrposc_1.0 <- (1 - bounds$pc_1.0) / (1 - bounds$pc_0.0)
bounds$lrposc_1.2 <- (1 - bounds$pc_1.2) / (1 - bounds$pc_0.0)
bounds$lrposc_1.4 <- (1 - bounds$pc_1.4) / (1 - bounds$pc_0.0)
bounds$lrposc_1.6 <- (1 - bounds$pc_1.6) / (1 - bounds$pc_0.0)
bounds$lrposc_1.8 <- (1 - bounds$pc_1.8) / (1 - bounds$pc_0.0)
bounds$lrposc_2.0 <- (1 - bounds$pc_2.0) / (1 - bounds$pc_0.0)
bounds$lrposc_2.2 <- (1 - bounds$pc_2.2) / (1 - bounds$pc_0.0)
bounds$lrposc_2.4 <- (1 - bounds$pc_2.4) / (1 - bounds$pc_0.0)
bounds$lrposc_2.6 <- (1 - bounds$pc_2.6) / (1 - bounds$pc_0.0)
bounds$lrposc_2.8 <- (1 - bounds$pc_2.8) / (1 - bounds$pc_0.0)
bounds$lrposc_3.0 <- (1 - bounds$pc_3.0) / (1 - bounds$pc_0.0)
bounds$lrnega_0.2 <- bounds$pa_0.2 / bounds$pa_0.0
bounds$lrnega_0.4 <- bounds$pa_0.4 / bounds$pa_0.0
bounds$lrnega_0.6 <- bounds$pa_0.6 / bounds$pa_0.0
bounds$lrnega_0.8 <- bounds$pa_0.8 / bounds$pa_0.0
bounds$lrnega_1.0 <- bounds$pa_1.0 / bounds$pa_0.0
bounds$lrnega_1.2 <- bounds$pa_1.2 / bounds$pa_0.0
bounds$lrnega_1.4 <- bounds$pa_1.4 / bounds$pa_0.0
bounds$lrnega_1.6 <- bounds$pa_1.6 / bounds$pa_0.0
bounds$lrnega_1.8 <- bounds$pa_1.8 / bounds$pa_0.0
bounds$lrnega_2.0 <- bounds$pa_2.0 / bounds$pa_0.0
bounds$lrnega_2.2 <- bounds$pa_2.2 / bounds$pa_0.0
bounds$lrnega_2.4 <- bounds$pa_2.4 / bounds$pa_0.0
bounds$lrnega_2.6 <- bounds$pa_2.6 / bounds$pa_0.0
bounds$lrnega_2.8 <- bounds$pa_2.8 / bounds$pa_0.0
bounds$lrnega_3.0 <- bounds$pa_3.0 / bounds$pa_0.0
bounds$lrnegb_0.2 <- bounds$pb_0.2 / bounds$pb_0.0
bounds$lrnegb_0.4 <- bounds$pb_0.4 / bounds$pb_0.0
bounds$lrnegb_0.6 <- bounds$pb_0.6 / bounds$pb_0.0
bounds$lrnegb_0.8 <- bounds$pb_0.8 / bounds$pb_0.0
bounds$lrnegb_1.0 <- bounds$pb_1.0 / bounds$pb_0.0
bounds$lrnegb_1.2 <- bounds$pb_1.2 / bounds$pb_0.0
bounds$lrnegb_1.4 <- bounds$pb_1.4 / bounds$pb_0.0
bounds$lrnegb_1.6 <- bounds$pb_1.6 / bounds$pb_0.0
bounds$lrnegb_1.8 <- bounds$pb_1.8 / bounds$pb_0.0
bounds$lrnegb_2.0 <- bounds$pb_2.0 / bounds$pb_0.0
bounds$lrnegb_2.2 <- bounds$pb_2.2 / bounds$pb_0.0
bounds$lrnegb_2.4 <- bounds$pb_2.4 / bounds$pb_0.0
bounds$lrnegb_2.6 <- bounds$pb_2.6 / bounds$pb_0.0
bounds$lrnegb_2.8 <- bounds$pb_2.8 / bounds$pb_0.0
bounds$lrnegb_3.0 <- bounds$pb_3.0 / bounds$pb_0.0
bounds$lrnegc_0.2 <- bounds$pc_0.2 / bounds$pc_0.0
bounds$lrnegc_0.4 <- bounds$pc_0.4 / bounds$pc_0.0
bounds$lrnegc_0.6 <- bounds$pc_0.6 / bounds$pc_0.0
bounds$lrnegc_0.8 <- bounds$pc_0.8 / bounds$pc_0.0
bounds$lrnegc_1.0 <- bounds$pc_1.0 / bounds$pc_0.0
bounds$lrnegc_1.2 <- bounds$pc_1.2 / bounds$pc_0.0
bounds$lrnegc_1.4 <- bounds$pc_1.4 / bounds$pc_0.0
bounds$lrnegc_1.6 <- bounds$pc_1.6 / bounds$pc_0.0
bounds$lrnegc_1.8 <- bounds$pc_1.8 / bounds$pc_0.0
bounds$lrnegc_2.0 <- bounds$pc_2.0 / bounds$pc_0.0
bounds$lrnegc_2.2 <- bounds$pc_2.2 / bounds$pc_0.0
bounds$lrnegc_2.4 <- bounds$pc_2.4 / bounds$pc_0.0
bounds$lrnegc_2.6 <- bounds$pc_2.6 / bounds$pc_0.0
bounds$lrnegc_2.8 <- bounds$pc_2.8 / bounds$pc_0.0
bounds$lrnegc_3.0 <- bounds$pc_3.0 / bounds$pc_0.0

# save all except large objects:
# save(list = ls()[substr(ls(), 1, 2) != "cr"], file = "crossrunbox1.Rdata")
