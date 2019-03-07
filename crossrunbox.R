library(Rmpfr)
library(crossrun)

# function for box with lowest probability for the target shift, among boxes
# with probability >=  target for shift 0. only boxes with positive corner
# probability considered. subsequent deletion within the border if possible is
# placed in separate function taken c and l for this box as input:
bestbox <- function(pt0    = crs100_0.0,
                    pts    = crs100_0.8,
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
  boxprt  <-
    two * (multm ^ (n1 - 1)) # initialize to impossible high value
  
  for (cc in 0:(n1 - 1))
    for (ll in 1:n1) {
      if (pt0n[cc + 1, ll] > nill &
          bpt0[cc + 1, ll] >= targt &
          bpttarg[cc + 1, ll] < boxprt) {
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
cutbox <- function(pt0    = crs100_0.0,
                   pts    = crs100_0.8,
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
  targt     <-
    targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n      <- pt0[[n1]]
  ptsn      <- pts[[n1]]
  bpt0      <-
    boxprobt(pt0n) # box probabilities for no shift, pt scale
  boxpt0    <-
    bpt0[c1 + 1, l1] # no shift probability of actual box, pt scale
  cornerpt0 <-
    pt0n[c1 + 1, l1] # no shift corner probability, pt scale
  finished  <- FALSE
  cbord     <- NA
  lbord     <- NA
  
  if (boxpt0 - cornerpt0 >= targt) {
    cutboxpt0 <-
      boxpt0 - cornerpt0 # pt of cutted box after removed corner
    cbord     <- c1 + 1
    lbord     <- l1 - 1
    
    while (!finished) {
      pt0n.directionc <- pt0n[cbord + 1, l1]
      pt0n.directionl <- pt0n[c1 + 1, lbord]
      ptsn.directionc <- ptsn[cbord + 1, l1]
      ptsn.directionl <- ptsn[c1 + 1, lbord]
      if ((cutboxpt0 - pt0n.directionc < targt |
           pt0n.directionc == 0) &
          (cutboxpt0 - pt0n.directionl < targt |
           pt0n.directionl == 0)) {
        finished <- TRUE
      } else if (cutboxpt0 - pt0n.directionc < targt |
                 pt0n.directionc == 0) {
        lstrip    <- pt0n[c1 + 1, lbord:1]
        nlstrip   <- length(lstrip)
        maxlstrip <- max((1:nlstrip)[lstrip > 0])
        lstrip    <- lstrip[(1:nlstrip) <= maxlstrip]
        lstripcum <- cumsum(lstrip)
        if (cutboxpt0 - max(lstripcum) >= targt) {
          lbord <- 0
        } else {
          # 0 cannot occurr
          lbord <-
            lbord + 1 - min((1:nlstrip)[cutboxpt0 - lstripcum < targt])
        }
        finished <- TRUE
      } else if (cutboxpt0 - pt0n.directionl < targt |
                 pt0n.directionl == 0) {
        cstrip    <- pt0n[(cbord + 1):n1, l1]
        ncstrip   <- length(cstrip)
        maxcstrip <- max((1:ncstrip)[cstrip > 0])
        cstrip    <- cstrip[(1:ncstrip) <= maxcstrip]
        cstripcum <- cumsum(cstrip)
        if (cutboxpt0 - max(cstripcum) >= targt) {
          cbord <- n1
        } else {
          # n1 cannot occurr
          cbord <-
            cbord + min((1:ncstrip)[cutboxpt0 - cstripcum < targt]) - 1
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
crs100_0.0 <- crossrunsymm(100, printn = TRUE)$pt

# shift 0.2 to 3:
crs100_0.2 <- crossrunshift(shift = 0.2, printn = TRUE)$pt
crs100_0.4 <- crossrunshift(shift = 0.4, printn = TRUE)$pt
crs100_0.6 <- crossrunshift(shift = 0.6, printn = TRUE)$pt
crs100_0.8 <- crossrunshift(shift = 0.8, printn = TRUE)$pt
crs100_1.0 <- crossrunshift(shift = 1.0, printn = TRUE)$pt
crs100_1.2 <- crossrunshift(shift = 1.2, printn = TRUE)$pt
crs100_1.4 <- crossrunshift(shift = 1.4, printn = TRUE)$pt
crs100_1.6 <- crossrunshift(shift = 1.6, printn = TRUE)$pt
crs100_1.8 <- crossrunshift(shift = 1.8, printn = TRUE)$pt
crs100_2.0 <- crossrunshift(shift = 2.0, printn = TRUE)$pt
crs100_2.2 <- crossrunshift(shift = 2.2, printn = TRUE)$pt
crs100_2.4 <- crossrunshift(shift = 2.4, printn = TRUE)$pt
crs100_2.6 <- crossrunshift(shift = 2.6, printn = TRUE)$pt
crs100_2.8 <- crossrunshift(shift = 2.8, printn = TRUE)$pt
crs100_3.0 <- crossrunshift(shift = 3.0, printn = TRUE)$pt

# Table 1 in "Run charts revisited", PLOS ONE November 25, 2014:
bounds <- data.frame(
  n = 10:100,
  ca = qbinom(0.05, 10:100 - 1, 0.5),
  la = round(log2(10:100) + 3)
)
row.names(bounds) <- bounds$n

# find alternative ("best") boxes:
bounds$cb <- NA
bounds$lb <- NA

for (nn in 10:100) {
  print(nn)
  bounds[bounds$n == nn, c("cb", "lb")] <- bestbox(n1 = nn)
}

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

# find no signal probabilities for the shifts, original and bestbox rules:
# initialization:
prec.use <- 120
one <- mpfr(1, prec.use)
two <- mpfr(2, prec.use)
mone <- mpfr(-1, prec.use)

# initialize to impossible (negative) value:
pat <- mpfr2array(rep(mone, 91 * 16), dim = c(91, 16))
rownames(pat) <- 10:100
colnames(pat) <-
  paste0("pat_", rep(c(0, 1, 2, 3), c(5, 5, 5, 1)), ".",
         c(rep(c(0, 2, 4, 6, 8), 3), 0))
pbt <- mpfr2array(rep(mone, 91 * 16), dim = c(91, 16))
rownames(pbt) <- 10:100
colnames(pbt) <- paste0("pbt_", rep(c(0, 1, 2, 3), c(5, 5, 5, 1)), ".",
                        c(rep(c(0, 2, 4, 6, 8), 3), 0))
pct <- mpfr2array(rep(mone, 91 * 16), dim = c(91, 16))
rownames(pct) <- 10:100
colnames(pct) <- paste0("pct_", rep(c(0, 1, 2, 3), c(5, 5, 5, 1)), ".",
                        c(rep(c(0, 2, 4, 6, 8), 3), 0))

# original boxes:
for (nn in 10:100) {
  print(nn)
  ca1 <- bounds$ca[bounds$n == nn]
  la1 <- bounds$la[bounds$n == nn]
  pat[nn - 9, 1] <- sum(crs100_0.0[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 2] <- sum(crs100_0.2[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 3] <- sum(crs100_0.4[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 4] <- sum(crs100_0.6[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 5] <- sum(crs100_0.8[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 6] <- sum(crs100_1.0[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 7] <- sum(crs100_1.2[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 8] <- sum(crs100_1.4[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 9] <- sum(crs100_1.6[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 10] <- sum(crs100_1.8[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 11] <- sum(crs100_2.0[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 12] <- sum(crs100_2.2[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 13] <- sum(crs100_2.4[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 14] <- sum(crs100_2.6[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 15] <- sum(crs100_2.8[[nn]][(ca1 + 1):nn, 1:la1])
  pat[nn - 9, 16] <- sum(crs100_3.0[[nn]][(ca1 + 1):nn, 1:la1])
}

# best boxes:
for (nn in 10:100) {
  print(nn)
  cb1 <- bounds$cb[bounds$n == nn]
  lb1 <- bounds$lb[bounds$n == nn]
  pbt[nn - 9, 1] <- sum(crs100_0.0[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 2] <- sum(crs100_0.2[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 3] <- sum(crs100_0.4[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 4] <- sum(crs100_0.6[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 5] <- sum(crs100_0.8[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 6] <- sum(crs100_1.0[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 7] <- sum(crs100_1.2[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 8] <- sum(crs100_1.4[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 9] <- sum(crs100_1.6[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 10] <- sum(crs100_1.8[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 11] <- sum(crs100_2.0[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 12] <- sum(crs100_2.2[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 13] <- sum(crs100_2.4[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 14] <- sum(crs100_2.6[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 15] <- sum(crs100_2.8[[nn]][(cb1 + 1):nn, 1:lb1])
  pbt[nn - 9, 16] <- sum(crs100_3.0[[nn]][(cb1 + 1):nn, 1:lb1])
}

# cut boxes:
pct <- pbt
colnames(pct) <- paste0("pct_", rep(c(0, 1, 2, 3), c(5, 5, 5, 1)), ".",
                        c(rep(c(0, 2, 4, 6, 8), 3), 0))

for (nn in 10:100) {
  print(nn)
  cb1    <- bounds$cb[bounds$n == nn]
  lb1    <- bounds$lb[bounds$n == nn]
  cbord1 <- bounds$cbord[bounds$n == nn]
  lbord1 <- bounds$lbord[bounds$n == nn]
  if (is.na(cbord1) == 0) {
    pct[nn - 9, 1] <- sum(crs100_0.0[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_0.0[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_0.0[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 2] <-
      sum(crs100_0.2[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_0.2[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_0.2[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 3] <-
      sum(crs100_0.4[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_0.4[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_0.4[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 4] <-
      sum(crs100_0.6[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_0.6[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_0.6[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 5] <-
      sum(crs100_0.8[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_0.8[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_0.8[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 6] <-
      sum(crs100_1.0[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_1.0[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_1.0[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 7] <-
      sum(crs100_1.2[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_1.2[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_1.2[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 8] <-
      sum(crs100_1.4[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_1.4[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_1.4[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 9] <-
      sum(crs100_1.6[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_1.6[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_1.6[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 10] <-
      sum(crs100_1.8[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_1.8[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_1.8[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 11] <-
      sum(crs100_2.0[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_2.0[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_2.0[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 12] <-
      sum(crs100_2.2[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_2.2[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_2.2[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 13] <-
      sum(crs100_2.4[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_2.4[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_2.4[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 14] <-
      sum(crs100_2.6[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_2.6[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_2.6[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 15] <-
      sum(crs100_2.8[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_2.8[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_2.8[[nn]][cb1 + 1, 1:lbord1])
    pct[nn - 9, 16] <-
      sum(crs100_3.0[[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
      sum(crs100_3.0[[nn]][(cbord1 + 1):nn, lb1]) +
      sum(crs100_3.0[[nn]][cb1 + 1, 1:lbord1])
  }
}

# computation of likelihood ratios:
# initialization to impossible (negative) value:
lrposa <- mpfr2array(rep(mone, 91 * 15), dim = c(91, 15))
rownames(lrposa) <- 10:100
colnames(lrposa) <-
  paste0("lrposa_", rep(c(0, 1, 2, 3), c(4, 5, 5, 1)), ".",
         rep(c(2, 4, 6, 8, 0), 3))
lrposb <- mpfr2array(rep(mone, 91 * 15), dim = c(91, 15))
rownames(lrposb) <- 10:100
colnames(lrposb) <-
  paste0("lrposb_", rep(c(0, 1, 2, 3), c(4, 5, 5, 1)), ".",
         rep(c(2, 4, 6, 8, 0), 3))
lrposc <- mpfr2array(rep(mone, 91 * 15), dim = c(91, 15))
rownames(lrposc) <- 10:100
colnames(lrposc) <-
  paste0("lrposc_", rep(c(0, 1, 2, 3), c(4, 5, 5, 1)), ".",
         rep(c(2, 4, 6, 8, 0), 3))
lrnega <- mpfr2array(rep(mone, 91 * 15), dim = c(91, 15))
rownames(lrnega) <- 10:100
colnames(lrnega) <-
  paste0("lrnega_", rep(c(0, 1, 2, 3), c(4, 5, 5, 1)), ".",
         rep(c(2, 4, 6, 8, 0), 3))
lrnegb <- mpfr2array(rep(mone, 91 * 15), dim = c(91, 15))
rownames(lrnegb) <- 10:100
colnames(lrnegb) <-
  paste0("lrnegb_", rep(c(0, 1, 2, 3), c(4, 5, 5, 1)), ".",
         rep(c(2, 4, 6, 8, 0), 3))
lrnegc <- mpfr2array(rep(mone, 91 * 15), dim = c(91, 15))
rownames(lrnegc) <- 10:100
colnames(lrnegc) <-
  paste0("lrnegc_", rep(c(0, 1, 2, 3), c(4, 5, 5, 1)), ".",
         rep(c(2, 4, 6, 8, 0), 3))

# computation:
lrposa[, "lrposa_0.2"] <-
  (two ^ (9:99) - pat[, "pat_0.2"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_0.4"] <-
  (two ^ (9:99) - pat[, "pat_0.4"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_0.6"] <-
  (two ^ (9:99) - pat[, "pat_0.6"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_0.8"] <-
  (two ^ (9:99) - pat[, "pat_0.8"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_1.0"] <-
  (two ^ (9:99) - pat[, "pat_1.0"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_1.2"] <-
  (two ^ (9:99) - pat[, "pat_1.2"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_1.4"] <-
  (two ^ (9:99) - pat[, "pat_1.4"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_1.6"] <-
  (two ^ (9:99) - pat[, "pat_1.6"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_1.8"] <-
  (two ^ (9:99) - pat[, "pat_1.8"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_2.0"] <-
  (two ^ (9:99) - pat[, "pat_2.0"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_2.2"] <-
  (two ^ (9:99) - pat[, "pat_2.2"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_2.4"] <-
  (two ^ (9:99) - pat[, "pat_2.4"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_2.6"] <-
  (two ^ (9:99) - pat[, "pat_2.6"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_2.8"] <-
  (two ^ (9:99) - pat[, "pat_2.8"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposa[, "lrposa_3.0"] <-
  (two ^ (9:99) - pat[, "pat_3.0"]) / (two ^ (9:99) - pat[, "pat_0.0"])
lrposb[, "lrposb_0.2"] <-
  (two ^ (9:99) - pbt[, "pbt_0.2"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_0.4"] <-
  (two ^ (9:99) - pbt[, "pbt_0.4"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_0.6"] <-
  (two ^ (9:99) - pbt[, "pbt_0.6"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_0.8"] <-
  (two ^ (9:99) - pbt[, "pbt_0.8"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_1.0"] <-
  (two ^ (9:99) - pbt[, "pbt_1.0"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_1.2"] <-
  (two ^ (9:99) - pbt[, "pbt_1.2"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_1.4"] <-
  (two ^ (9:99) - pbt[, "pbt_1.4"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_1.6"] <-
  (two ^ (9:99) - pbt[, "pbt_1.6"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_1.8"] <-
  (two ^ (9:99) - pbt[, "pbt_1.8"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_2.0"] <-
  (two ^ (9:99) - pbt[, "pbt_2.0"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_2.2"] <-
  (two ^ (9:99) - pbt[, "pbt_2.2"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_2.4"] <-
  (two ^ (9:99) - pbt[, "pbt_2.4"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_2.6"] <-
  (two ^ (9:99) - pbt[, "pbt_2.6"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_2.8"] <-
  (two ^ (9:99) - pbt[, "pbt_2.8"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposb[, "lrposb_3.0"] <-
  (two ^ (9:99) - pbt[, "pbt_3.0"]) / (two ^ (9:99) - pbt[, "pbt_0.0"])
lrposc[, "lrposc_0.2"] <-
  (two ^ (9:99) - pct[, "pct_0.2"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_0.4"] <-
  (two ^ (9:99) - pct[, "pct_0.4"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_0.6"] <-
  (two ^ (9:99) - pct[, "pct_0.6"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_0.8"] <-
  (two ^ (9:99) - pct[, "pct_0.8"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_1.0"] <-
  (two ^ (9:99) - pct[, "pct_1.0"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_1.2"] <-
  (two ^ (9:99) - pct[, "pct_1.2"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_1.4"] <-
  (two ^ (9:99) - pct[, "pct_1.4"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_1.6"] <-
  (two ^ (9:99) - pct[, "pct_1.6"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_1.8"] <-
  (two ^ (9:99) - pct[, "pct_1.8"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_2.0"] <-
  (two ^ (9:99) - pct[, "pct_2.0"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_2.2"] <-
  (two ^ (9:99) - pct[, "pct_2.2"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_2.4"] <-
  (two ^ (9:99) - pct[, "pct_2.4"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_2.6"] <-
  (two ^ (9:99) - pct[, "pct_2.6"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_2.8"] <-
  (two ^ (9:99) - pct[, "pct_2.8"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrposc[, "lrposc_3.0"] <-
  (two ^ (9:99) - pct[, "pct_3.0"]) / (two ^ (9:99) - pct[, "pct_0.0"])
lrnega[, "lrnega_0.2"] <- pat[, "pat_0.2"] / pat[, "pat_0.0"]
lrnega[, "lrnega_0.4"] <- pat[, "pat_0.4"] / pat[, "pat_0.0"]
lrnega[, "lrnega_0.6"] <- pat[, "pat_0.6"] / pat[, "pat_0.0"]
lrnega[, "lrnega_0.8"] <- pat[, "pat_0.8"] / pat[, "pat_0.0"]
lrnega[, "lrnega_1.0"] <- pat[, "pat_1.0"] / pat[, "pat_0.0"]
lrnega[, "lrnega_1.2"] <- pat[, "pat_1.2"] / pat[, "pat_0.0"]
lrnega[, "lrnega_1.4"] <- pat[, "pat_1.4"] / pat[, "pat_0.0"]
lrnega[, "lrnega_1.6"] <- pat[, "pat_1.6"] / pat[, "pat_0.0"]
lrnega[, "lrnega_1.8"] <- pat[, "pat_1.8"] / pat[, "pat_0.0"]
lrnega[, "lrnega_2.0"] <- pat[, "pat_2.0"] / pat[, "pat_0.0"]
lrnega[, "lrnega_2.2"] <- pat[, "pat_2.2"] / pat[, "pat_0.0"]
lrnega[, "lrnega_2.4"] <- pat[, "pat_2.4"] / pat[, "pat_0.0"]
lrnega[, "lrnega_2.6"] <- pat[, "pat_2.6"] / pat[, "pat_0.0"]
lrnega[, "lrnega_2.8"] <- pat[, "pat_2.8"] / pat[, "pat_0.0"]
lrnega[, "lrnega_3.0"] <- pat[, "pat_3.0"] / pat[, "pat_0.0"]
lrnegb[, "lrnegb_0.2"] <- pbt[, "pbt_0.2"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_0.4"] <- pbt[, "pbt_0.4"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_0.6"] <- pbt[, "pbt_0.6"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_0.8"] <- pbt[, "pbt_0.8"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_1.0"] <- pbt[, "pbt_1.0"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_1.2"] <- pbt[, "pbt_1.2"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_1.4"] <- pbt[, "pbt_1.4"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_1.6"] <- pbt[, "pbt_1.6"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_1.8"] <- pbt[, "pbt_1.8"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_2.0"] <- pbt[, "pbt_2.0"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_2.2"] <- pbt[, "pbt_2.2"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_2.4"] <- pbt[, "pbt_2.4"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_2.6"] <- pbt[, "pbt_2.6"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_2.8"] <- pbt[, "pbt_2.8"] / pbt[, "pbt_0.0"]
lrnegb[, "lrnegb_3.0"] <- pbt[, "pbt_3.0"] / pbt[, "pbt_0.0"]
lrnegc[, "lrnegc_0.2"] <- pct[, "pct_0.2"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_0.4"] <- pct[, "pct_0.4"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_0.6"] <- pct[, "pct_0.6"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_0.8"] <- pct[, "pct_0.8"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_1.0"] <- pct[, "pct_1.0"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_1.2"] <- pct[, "pct_1.2"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_1.4"] <- pct[, "pct_1.4"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_1.6"] <- pct[, "pct_1.6"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_1.8"] <- pct[, "pct_1.8"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_2.0"] <- pct[, "pct_2.0"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_2.2"] <- pct[, "pct_2.2"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_2.4"] <- pct[, "pct_2.4"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_2.6"] <- pct[, "pct_2.6"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_2.8"] <- pct[, "pct_2.8"] / pct[, "pct_0.0"]
lrnegc[, "lrnegc_3.0"] <- pct[, "pct_3.0"] / pct[, "pct_0.0"]

# include probability information in extended bounds:
pa <- pat
pb <- pbt
pc <- pct
prec.use <- 120
one <- mpfr(1, prec.use)
two <- mpfr(2, prec.use)
mone <- mpfr(-1, prec.use)

for (var in 1:16) pa[,var] <- pat[,var]/(two^(9:99))
for (var in 1:16) pb[,var] <- pbt[,var]/(two^(9:99))
for (var in 1:16) pc[,var] <- pct[,var]/(two^(9:99))

boundspl <- cbind(
  bounds, 
  asNumeric(pa), asNumeric(pb), asNumeric(pc),
  asNumeric(lrposa), asNumeric(lrposb), asNumeric(lrposc),
  asNumeric(lrnega), asNumeric(lrnegb), asNumeric(lrnegc)) 

names(boundspl)[8:55] <- paste0(
  substr(names(boundspl)[8:55],1,2),
  substr(names(boundspl)[8:55],4,7))

# save objects:
save(bounds, boundspl,
     pat, pbt, pct,
     pa, pb, pc,
     lrposa, lrposb, lrposc,
     lrnega, lrnegb, lrnegc,
     file = "data/boundspl.Rdata")

save(list = ls(pattern = 'crs100_'), 
     file="data/crs100.RData")

# Save data for box figures
x <- list(
  x0 = lapply(crs100_0.0, asNumeric),
  x1 = lapply(crs100_1.0, asNumeric),
  x2 = lapply(crs100_2.0, asNumeric))

saveRDS(x, 'data/x.rds')

# save(lrposa, lrposb, lrposc, lrnega, lrnegb, lrnegc, file = "lr.RData")
# save(bounds, file = "bounds.Rdata")
# save(pat, pbt, pct, file = "pt.RData")
