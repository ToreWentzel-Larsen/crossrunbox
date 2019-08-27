# ###############################################################################
# This script reproduces the data objects and figures from the Smooth Operator
# article by Jacob Anhøj and Tore Wentzel-Larsen, R Journal 2019.
#
# The objects of interest are cr_dists and cr_bounds.
#
# cr_dists: a list with probabilities for the joint distribution of number of
# crossings (C) and longest run (L) in multiple precision (mpfr) format. To get
# the matrix of probabilities for, say, N = 11 and no shift (SD = 0), use
# cr_dists$pt_0.0[[11]].
#
# cr_bounds: a data frame with critical values for longest run and number of
# crossings together with probabilities and log-likelihood ratios for the Anhøj,
# best box and cut box rules.
# Variables:
#   ca = lower limit for number of crossings, Anhøj rules. 
#   la = upper limit for longest run, Anhøj rules. 
#   cb = lower limit for number of crossings, best box rules. 
#   lb = upper limit for longest run, best box rules. 
#   cbord/lbord = coordinates for the cut box adjustment.
#   pa_n.m/pb_n.m / pcbord_n.m = probality of no signal. n.m = size of the shift
#                              in standard deviation units.
#   loglrpos_n.m / loglrneg_n.m = positive and negative log-likehood ratios.
# 
# For the sake of speed and memory consumption the scrip is by default set to
# produce output for N = 10-40 and SD = 0-2. To reproduce all data from the
# article, change the parameters nmax and smax to 100 and 3 respectively.
# 
# Jacob Anhøj & Tore Wentzel-Larsen 24 Mar 2019
################################################################################

# Load libraries ----
library(Rmpfr)
library(crossrun)
library(tidyverse)

# Set parameters ----
#nmax         <- 40      # Max N to include in computations.
# crossrunreg: changed to 100
nmax         <- 100
smax         <- 2       # Max shift in SD units to include in computations.
target       <- 0.925   # Target specificity for best box and cut box.
target_shift <- 0.8     # Target shift for best box and cut box.

## Probably no need to change anything below this line
nmin     <- 10
shifts   <- seq(0, smax, by = 0.2)

# bestbox function ----
## Function for box with lowest probability for the target shift, among boxes
## with probability >= target for shift = 0.
bestbox <- function(pt0,
                    pts,
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
  bpt0    <- boxprobt(pt0n)  # box probabilities for no shift
  bpttarg <- boxprobt(ptsn)  # box probabilities for target shift
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
    }
  return(c(c1, l1))
} # end function bestbox

# cutbox function ----
## Function for cutting a box while keeping probability >= target for shift = 0.
## No cutting if the corner cannot be removed.
cutbox <- function(pt0,
                   pts,
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
  bpt0      <-
    boxprobt(pt0n)   # box probabilities for no shift, pt scale
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
    }
  }
  return(c(cbord, lbord))
} # end function cutbox

# crs function to compute joint distributions of C and L ----
crs <- function(nmax = 12,
                shifts = seq(0, 2, by = 0.2)) {
  
  crs <- list()
  for (s in shifts) {
    r <- paste0('pt_', format(s, nsmall = 1))
    print(paste('Joint distribution:', r))
    crs[[r]] <- crossrunshift(nmax, s)$pt
  }
  return(crs)
} # End crs function

# bounds function to compute limits and diagnostics for runs rules ----
bounds <- function(crs, target = 0.925, target_shift = 0.8) {
  nmin     <- 10
  nmax     <- length(crs[[1]])
  shiftsc  <- format(shifts, nsmall = 1)
  nshifts  <- length(shifts)
  prec.use <- 120
  one      <- mpfr(1, prec.use)
  two      <- mpfr(2, prec.use)
  mone     <- mpfr(-1, prec.use)
  pt0      <- crs$pt_0.0
  # pts      <- crs$pt_0.8
  pts      <- crs[[paste0('pt_', format(target_shift, nsmall = 1))]]
  
  # Begin bounds table
  bounds <- data.frame(
    n  = nmin:nmax,
    ca = qbinom(0.05, nmin:nmax - 1, 0.5),
    la = round(log2(nmin:nmax) + 3)
  )
  row.names(bounds) <- bounds$n
  
  # find best boxes
  bounds$cb <- NA
  bounds$lb <- NA
  
  for (nn in nmin:nmax) {
    print(paste('bestbox:', nn))
    bounds[bounds$n == nn, c('cb', 'lb')] <- bestbox(pt0,
                                                     pts,
                                                     n1 = nn,
                                                     target = target)
  }
  
  # find cut  boxes
  bounds$cbord <- NA
  bounds$lbord <- NA
  
  for (nn in nmin:nmax) {
    print(paste('cutbox', nn))
    bounds[bounds$n == nn, c('cbord', 'lbord')] <-
      cutbox(pt0,
             pts,
             n1 = nn,
             target = target,
             c1 = bounds$cb[bounds$n == nn],
             l1 = bounds$lb[bounds$n == nn])
  }
  
  # Find no signal probabilities
  ## initialize to impossible (negative) value:
  pat <- mpfr2array(rep(mone, (nmax - 9) * nshifts),
                    dim = c((nmax - 9), nshifts),
                    dimnames = list(nmin:nmax, NULL))
  
  pbt       <- pat
  pct       <- pat
  loglrposa <- pat[ , -1]
  loglrposb <- loglrposa
  loglrposc <- loglrposa
  loglrnega <- loglrposa
  loglrnegb <- loglrposa
  loglrnegc <- loglrposa
  
  colnames(pat)       <- paste0('pat_', shiftsc)
  colnames(pbt)       <- paste0('pbt_', shiftsc)
  colnames(pct)       <- paste0('pct_', shiftsc)
  colnames(loglrposa) <- paste0('loglrposa_', shiftsc[-1])
  colnames(loglrposb) <- paste0('loglrposb_', shiftsc[-1])
  colnames(loglrposc) <- paste0('loglrposc_', shiftsc[-1])
  colnames(loglrnega) <- paste0('loglrnega_', shiftsc[-1])
  colnames(loglrnegb) <- paste0('loglrnegb_', shiftsc[-1])
  colnames(loglrnegc) <- paste0('loglrnegc_', shiftsc[-1])
  
  ## calculations
  for (nn in nmin:nmax) {
    ca1    <- bounds$ca[bounds$n == nn]
    la1    <- bounds$la[bounds$n == nn]
    cb1    <- bounds$cb[bounds$n == nn]
    lb1    <- bounds$lb[bounds$n == nn]
    cbord1 <- bounds$cbord[bounds$n == nn]
    lbord1 <- bounds$lbord[bounds$n == nn]
    
    for (s in shifts) {
      i              <- match(s, shifts)
      p              <- format(s, nsmall = 1)
      p              <- paste0('pt_', p)
      pat[nn - 9, i] <- sum(crs[[p]][[nn]][(ca1 + 1):nn, 1:la1])
      pbt[nn - 9, i] <- sum(crs[[p]][[nn]][(cb1 + 1):nn, 1:lb1])
      pct[nn - 9, i] <- pbt[nn - 9, i]
      
      if (!is.na(cbord1)) {
        pct[nn - 9, i] <- 
          sum(crs[[p]][[nn]][(cb1 + 2):nn, 1:(lb1 - 1)]) +
          sum(crs[[p]][[nn]][(cbord1 + 1):nn, lb1]) +
          sum(crs[[p]][[nn]][cb1 + 1, 1:lbord1])
      }
    }
  }
  
  # Find likelihood ratios
  for (s in shiftsc[shifts > 0]) {
    pats <- paste0('pat_', s)
    pbts <- paste0('pbt_', s)
    pcts <- paste0('pct_', s)
    loglrposa1 <- paste0('loglrposa_', s)
    loglrposb1 <- paste0('loglrposb_', s)
    loglrposc1 <- paste0('loglrposc_', s)
    loglrnega1 <- paste0('loglrnega_', s)
    loglrnegb1 <- paste0('loglrnegb_', s)
    loglrnegc1 <- paste0('loglrnegc_', s)
    
    loglrposa[, loglrposa1] <-
      log(two ^ (nmin:nmax - 1) - pat[, pats]) - 
      log(two ^ (nmin:nmax - 1) - pat[, 'pat_0.0'])
    loglrposb[, loglrposb1] <-
      log(two ^ (nmin:nmax - 1) - pbt[, pbts]) - 
      log(two ^ (nmin:nmax - 1) - pbt[, 'pbt_0.0'])
    loglrposc[, loglrposc1] <-
      log(two ^ (nmin:nmax - 1) - pct[, pcts]) - 
      log(two ^ (nmin:nmax - 1) - pct[, 'pct_0.0'])
    
    loglrnega[, loglrnega1] <-
      log(pat[, pats]) - log(pat[, 'pat_0.0'])
    loglrnegb[, loglrnegb1] <-
      log(pbt[, pbts]) - log(pbt[, 'pbt_0.0'])
    loglrnegc[, loglrnegc1] <-
      log(pct[, pcts]) - log(pct[, 'pct_0.0'])
  }
  
  # Finish bounds table including probability information
  pa <- pat / (two ^ (nmin:nmax - 1))
  pb <- pbt / (two ^ (nmin:nmax - 1))
  pc <- pct / (two ^ (nmin:nmax - 1))
  
  bounds <- cbind(
    bounds,
    asNumeric(pa),        asNumeric(pb),        asNumeric(pc),
    asNumeric(loglrposa), asNumeric(loglrposb), asNumeric(loglrposc),
    asNumeric(loglrnega), asNumeric(loglrnegb), asNumeric(loglrnegc)
  )
  
  ## Fix column names, pat -> pa etc.
  names(bounds) <- sub("(^p.{1}).", "\\1\\", names(bounds))
  
  return(bounds)
} # End bounds function

# crplot function to plot joint CL probabilites and box bounds ----
crplot <- function(bounds, cr_dists, n = 11, labels = T) {
  ca    <- bounds$ca[bounds$n == n]
  la    <- bounds$la[bounds$n == n]
  pa    <- bounds$pa_0.0[bounds$n == n]
  cb    <- bounds$cb[bounds$n == n]
  lb    <- bounds$lb[bounds$n == n]
  pb    <- bounds$pb_0.0[bounds$n == n]
  cbord <- bounds$cbord[bounds$n == n]
  lbord <- bounds$lbord[bounds$n == n]
  pc    <- bounds$pc_0.0[bounds$n == n]
  pt    <- paste0('pt', n)
  # cr    <- crs
  
  d <- cr_dists[[pt]] %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column('C') %>% 
    gather('L', 'times', -C) %>% 
    mutate(L = as.integer(L),
           C = as.integer(C),
           p = times / sum(times),
           y = times,
           col = (times / max(times)) < 0.5 & (times / max(times)) > 0)
  
  p <- ggplot(d, aes(L, C, 
                     fill = times > 0,
                     alpha = times / max(times))) +
    geom_raster() +
    geom_rect(aes(xmin = 0.5,          # Anhoej box
                  xmax = la + 0.5,
                  ymin = ca - 0.5,
                  ymax = max(C) + 0.5),
              colour = '#F8766D',
              size = 1,
              fill   = NA) +
    geom_rect(aes(xmin = 0.5,          # Best box
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = max(C) + 0.5),
              linetype = 2,
              size = 1,
              colour   = '#00BA38',
              fill = NA) +
    geom_rect(aes(xmin = lbord + 0.5,  # Cut box, horizontal
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = cb + 0.5),
              colour   = '#619CFF',
              fill     = NA,
              linetype = 4,
              size = 1,
              na.rm    = T) +
    geom_rect(aes(xmin = lb - 0.5,     # Cut box, vertical
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = cbord - 0.5), 
              colour   = '#619CFF',
              linetype = 4,
              fill     = NA,
              na.rm    = T)
  
  if(labels) {
    p <- p +
      geom_text(aes(label = y,
                    colour = col),
                alpha = 1,
                size = 3)
  }
  
  p <- p +
    scale_y_reverse(breaks = 0:max(d$C)) +
    scale_x_continuous(position = 'top',
                       breaks = 1:max(d$L)) +
    scale_fill_manual(values = c('white', 'black')) +
    scale_colour_manual(values = c('white', 'black')) +
    theme_minimal() +
    theme(axis.title.y = element_text(angle = 0, hjust = 0),
          axis.title.x = element_text(hjust = 0),
          panel.grid.major = element_blank(),
          aspect.ratio = 1,
          legend.position = 'none') +
    labs(caption = paste0('N =', n, '\n',
                          'Rules specificity: ',
                          'Anhøj = ', round(pa, 3), ', ',
                          'Best Box = ', round(pb, 3), ', ',
                          'Cut Box = ', round(pc, 3)))
  plot(p)
} # End crplot function

# Create data objects ----
## Joint distribution matrices
cr_dists    <- crs(nmax, shifts)

## Data frame with limits and diagnostics for runs rules
cr_bounds <- bounds(cr_dists, target, target_shift)

## Bounds data in tall format
cr_bounds_tall <- cr_bounds %>% 
  select(-(ca:lbord)) %>%
  gather('key', 'val', -n) %>% 
  separate(key, c('test', 'shift'), '_') %>% 
  mutate(rule = substring(test, nchar(test)),
         test = substring(test, 1, nchar(test) - 1),
         shift = as.numeric(shift)) %>% 
  mutate(rule = fct_recode(rule, 
                           anhoej = 'a',
                           `best box` = 'b',
                           `cut box` = 'c')) %>% 
  spread(test, val)

# Figures ----
## Plot joint distribution matrix
crplot(cr_bounds, map(cr_dists$pt_0.0, asNumeric), 11, labels = T)
crplot(cr_bounds, map(cr_dists$pt_0.0, asNumeric), 19, labels = T)

## Plot power function
ggplot(cr_bounds_tall, aes(n, 1 - p, colour = rule)) +
  geom_line(size = 1) +
  facet_wrap(~ shift, ncol = 4) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Power function',
       y = 'Probability of signal',
       x = 'N')

## Plot specificity
ggplot(filter(cr_bounds_tall, shift == 0), aes(n, p, colour = rule)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Specificity',
       y = 'Probability of true negative',
       x = 'N')

## Plot sensitivity for shift = 0.8 SD
ggplot(filter(cr_bounds_tall, shift == 0.8), aes(n, 1 - p, colour = rule)) +
  geom_line(size = 1) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Sensitivity (shift = 0.8 SD)',
       y = 'Probability of true positive',
       x = 'N')

## Plot LR+
ggplot(filter(cr_bounds_tall, !is.na(loglrpos)), 
       aes(n, exp(loglrpos), colour = rule)) +
  geom_line(size = 0.75) +
  geom_hline(yintercept = 10) +
  facet_wrap(~ shift, ncol = 5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Positive likelihood ratios',
       y = 'LR+',
       x = 'N')

## Plot LR-
ggplot(filter(cr_bounds_tall, !is.na(loglrpos)),
       aes(n, exp(loglrneg), colour = rule)) +
  geom_line(size = 0.75) +
  geom_hline(yintercept = 0.1) +
  facet_wrap(~ shift, ncol = 5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Negative likelihood ratios',
       y = 'LR-',
       x = 'N')


# tries building domains for declaring no special cause variation
# from lower left corner and up , for one n
domainbuild <- function(pt0,
                    pts,
                    target = 0.925,
                    n1     = 100,
                    mult   = 2,
                    prec   = 120,
                    check=0) {
  nill    <- mpfr(0, prec)
  one     <- mpfr(1, prec)
  two     <- mpfr(2, prec)
  multm   <- mpfr(mult, prec)
  targetm <- mpfr(target, prec)
  targt   <- targetm * (multm ^ (n1 - 1)) # target on "times" scale
  pt0n    <- pt0[[n1]]
  ptsn    <- pts[[n1]]
  bpt0    <- boxprobt(pt0n)  # box probabilities for no shift
  bpttarg <- boxprobt(ptsn)  # box probabilities for target shift
  # initializing, domain lower right cell, and
  # its (empty) vertical and horizontal strips
  domind <- Rmpfr::mpfr2array(rep(nill,n1*n1), dim = c(n1, n1))
  rownames(domind) <- 0:(n1-1)
  colnames(domind) <- 1:n1
  domind[,1] <- one
  domind[n1,] <- one
  pt0dom <- sum(domind*pt0n)
  ptsdom <- sum(domind*ptsn)
  step <- 0
  while(pt0dom<targt) {
    step <- step + 1
    domindold <- domind
    candid <- Rmpfr::mpfr2array(rep(nill,n1*n1), dim = c(n1, n1))
    for (cp in 1:(n1-1)) for (ll in 2:n1) 
      if (domindold[cp,ll]==0&domindold[cp+1,ll]==1&domindold[cp,ll-1]==1)
        candid[cp,ll] <- one
    ndom <- asNumeric(sum(candid))
    candidc <- rep(0,ndom)
    candidl <- rep(0,ndom)
    candidpt0 <- rep(nill, ndom)
    candidpts <- rep(nill, ndom)
    domi <- 1
    for (cc in 0:(n1-1)) for (ll in 1:n1) if (candid[cc+1,ll]==1) {
      candidc[domi] <- cc
      candidl[domi] <- ll
      candidpt0[domi] <- pt0n[cc+1,ll]
      candidpts[domi] <- ptsn[cc+1,ll]
      domi <- domi + 1
    }
    candidc <- candidc[order(candidpts)]
    candidl <- candidl[order(candidpts)]
    cnew <- candidc[1]
    lnew <- candidl[1]
    domind[cnew+1,lnew] <- one
    if (check==1) domindold <- domind
    if (max(pt0n[cnew+1,lnew:n1])==nill&domind[cnew+2,n1]==one)
      domind[cnew+1,lnew:n1] <- one
    if (max(pt0n[1:(cnew+1),lnew])==nill&domind[1,lnew-1]==one)
      domind[1:(cnew+1),lnew] <- one
    if (check==1) {
      print(step)
      domindeq <- max(abs(domindold-domind))==nill
      if (domindeq) print(asNumeric(domind))
      if (!domindeq) {
        domindcomb <- cbind(domindold,domind)
        rownames(domindcomb) <- 0:(n1-1)
        colnames(domindcomb) <- rep(1:n1,2)
        print(asNumeric(domindcomb))
        }
    } # end check output 
    pt0dom <- sum(domind*pt0n)
    ptsdom <- sum(domind*ptsn)
  } # end while below target
  for (cp in 1:n1) {
    rgcp <- domind[cp,]
    pt0ncp <- pt0n[cp,]
    if (max(pt0ncp[rgcp==one])==nill) 
      domind[cp,] <- nill
  } # end delete redundant rows
  for (ll in 1:n1) {
    rgll <- domind[,ll]
    pt0nll <- pt0n[,ll]
    if (max(pt0nll[rgll==one])==nill) 
      domind[,ll] <- nill
  } # end delete redundant columns
  return(domind)
} # end function domainbuild

# check in two cases:
# sequence length 11
d11 <- domainbuild(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                          n1=11, check=1) # 36 steps
asNumeric(d11) # sensitivity and specificity of rule d11:
sum(cr_dists[[1]][[11]]*d11)/sum(cr_dists[[1]][[11]])
1 - sum(cr_dists[[5]][[11]]*d11)/sum(cr_dists[[5]][[11]])
# sequence length 19
d19 <- domainbuild(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                n1=19, check=1) #  88 steps
asNumeric(d19) # sensitivity and specificity of rule d19:
sum(cr_dists[[1]][[19]]*d19)/sum(cr_dists[[1]][[19]])
1 - sum(cr_dists[[5]][[19]]*d19)/sum(cr_dists[[5]][[19]])
# details on candidates at step 68:
cr_dists[[5]][[19]][6+1,5]/sum(cr_dists[[5]][[19]])
cr_dists[[5]][[19]][8+1,6]/sum(cr_dists[[5]][[19]])
# probabilities in the symmetric case:
asNumeric(cr_dists[[1]][[19]])
# sequence length 25
d25 <- domainbuild(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                n1=25, check=1) # 135 steps
asNumeric(d25) # sensitivity and specificity of rule d25:
sum(cr_dists[[1]][[25]]*d25)/sum(cr_dists[[1]][[25]])
1 - sum(cr_dists[[5]][[25]]*d25)/sum(cr_dists[[5]][[25]])

# domain rules, function:
domainrulesf <- function(pt0,
                     pts,
                     target = 0.925,
                     nmin     = 10,
                     nmax     = 100,
                     mult   = 2,
                     prec   = 120,
                     printn=0) {
  one     <- mpfr(1, prec)
  if (printn==1) print(nmin)
  if (printn==1) print(Sys.time())
  dommatr <- list(domm=Rmpfr::mpfr2array(one, dim = c(nmin, nmin)))
  dommatr[[1]] <- domainbuild(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                           n1=nmin)
  for (nn in (nmin+1):nmax) {
    if (printn==1) print(nn)
    if (printn==1) print(Sys.time())
    dommatr[[nn+1-nmin]] <- domainbuild(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                                     n1=nn)
  }
  names(dommatr) <- nmin:nmax
  return(dommatr)
} # end function domainrulesf
  
domrules.10.50 <- domainrulesf(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                            nmax=50, printn=1)
domrules.51.60 <- domainrulesf(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                            nmin=51, nmax=60, printn=1)
# merge:
for (nn in 1:10) domrules.10.50[[41+nn]] <- domrules.51.60[[nn]]
names(domrules.10.50)[42:51] <- 51:60
# check:
sum(abs(domrules.10.50[[42]] - domrules.51.60[[1]]))
sum(abs(domrules.10.50[[43]] - domrules.51.60[[2]]))
sum(abs(domrules.10.50[[44]] - domrules.51.60[[3]]))
sum(abs(domrules.10.50[[45]] - domrules.51.60[[4]]))
sum(abs(domrules.10.50[[46]] - domrules.51.60[[5]]))
sum(abs(domrules.10.50[[47]] - domrules.51.60[[6]]))
sum(abs(domrules.10.50[[48]] - domrules.51.60[[7]]))
sum(abs(domrules.10.50[[49]] - domrules.51.60[[8]]))
sum(abs(domrules.10.50[[50]] - domrules.51.60[[9]]))
sum(abs(domrules.10.50[[51]] - domrules.51.60[[10]]))
# all 0, ok, delete domrules.51.60
rm(domrules.51.60)

domrules.61.70 <- domainrulesf(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                               nmin=61, nmax=70, printn=1)
save.image("M:/Statistisk prosesskontroll/crossrunboxarticle-master/crossrundomain.RData")

# merge:
domrules.10.70 <- domrules.10.50
for (nn in 1:10) domrules.10.70[[51+nn]] <- domrules.61.70[[nn]]
names(domrules.10.70)[52:61] <- 61:70
# check:
for (nn in 1:51) print(sum(abs(domrules.10.70[[nn]] - domrules.10.50[[nn]])))
for (n1 in 1:10) print(sum(abs(domrules.10.80[[61+n1]] - domrules.71.80[[n1]])))
rm(domrules.10.70,domrules.71.80)

domrules.71.80 <- domainrulesf(pt0=cr_dists[[1]], pts=cr_dists[[5]], 
                               nmin=71, nmax=80, printn=1)
save.image("M:/Statistisk prosesskontroll/crossrunboxarticle-master/crossrundomain.RData")

# merge:
domrules.10.80 <- domrules.10.70
for (nn in 1:10) domrules.10.80[[61+nn]] <- domrules.71.80[[nn]]
names(domrules.10.80)[62:71] <- 71:80
# check:
for (nn in 1:61) print(sum(abs(domrules.10.80[[nn]] - domrules.10.70[[nn]])))
for (n1 in 1:10) print(sum(abs(domrules.10.70[[51+n1]] - domrules.61.70[[n1]])))
rm(domrules.10.50,domrules.51.60,domrules.61.70)

# light version:
domruleslight.10.80 <- domrules.10.80
for (nn in 10:80) domruleslight.10.80[[nn-9]] <- asNumeric(domrules.10.80[[nn-9]])
# checks:
ndraw <- sample(10:80,1)
domruleslight.10.80[[ndraw-9]]

# simpler description of the rules
# auxiliary functioN:
lmaxrowf <- function(drulr) {
  max1 <- max(drulr)
  if (max1==0) res <- 0
  if (max1==1) {
    l1 <- 1:length(drulr)
    res <- max(l1[drulr==1])
  } # # end if l>0 occurs
  return(res)
} # end function lmaxrowf
# describe the rules by cmin, cmax and lmax for all occurring lmax:
cmicmalmaf <- function(drul, check=0) {
  n1 <- dim(drul)[1]
  c1 <- c(0:(n1-1))
  lmax1 <- apply(drul,1,lmaxrowf)
  if (check==1) print(lmax1)
  clm <- data.frame(c1=c1,lm=lmax1)
  clm <- clm[clm$lm>0,]
  if (check==1) print(clm)
  lmu <- sort(unique(clm$lm))
  if (check==1) print(lmu)
  lmun <- length(lmu)
  if (check==1) print(lmun)
  cmin1 <- rep(0,lmun)
  cmax1 <- rep(0,lmun)
  for (ll in c(1:lmun)) {
    cmin1[ll] <- min(clm$c1[clm$lm==lmu[ll]])
    cmax1[ll] <- max(clm$c1[clm$lm==lmu[ll]])
  } # end find cmin and cmax
  return(data.frame(cmin=cmin1,cmax=cmax1,l=lmu))
} # end function cmicmalmaf

# check:
ndraw <- sample(10:80,1)
domruleslight.10.80[[ndraw-9]]
cmicmalmaf(domruleslight.10.80[[ndraw-9]])

dommatr <- list(domm=Rmpfr::mpfr2array(one, dim = c(nmin, nmin)))


# create a list with the description of the rules:
domrulesdesc.10.80 <- list(cmicmalmaf(domruleslight.10.80[[1]]))
for (nn in 11:80) domrulesdesc.10.80[[nn-9]] <- 
  cmicmalmaf(domruleslight.10.80[[nn-9]])
names(domrulesdesc.10.80) <- 10:80

# stores information on the rules, without Mpfr precision:
save(domruleslight.10.80, domrulesdesc.10.80,
     file="domain build rules 10 to 80.Rdata")


# show rules:
for (nn in 10:80) {
  print(nn)
  print(domrulesdesc.10.80[[nn-9]])
}

# specificities:
cr_bounds$pd_0.0 <- rep(NA,91)
for (nn in 10:80) cr_bounds$pd_0.0[nn-9] <- 
  asNumeric(sum(cr_dists[[1]][[nn]]*domrules.10.80[[nn-9]])/sum(cr_dists[[1]][[nn]]))
# 1-sensitivities (shift 0.8):
cr_bounds$pd_0.8 <- rep(NA,91)
for (nn in 10:80) cr_bounds$pd_0.8[nn-9] <- 
  asNumeric(sum(cr_dists[[5]][[nn]]*domrules.10.80[[nn-9]])/sum(cr_dists[[5]][[nn]]))
1 - cr_bounds$pd_0.8
cr_bounds[1:71,c("pd_0.0","pd_0.8")]
# specificities and sensitivities:
spsecd <-  data.frame(n=cr_bounds$n[1:71],
           spes_c=cr_bounds$pc_0.0[1:71], spes_d=cr_bounds$pd_0.0[1:71],
           sens_c=1-cr_bounds$pc_0.8[1:71],
           sens_d=1-cr_bounds$pd_0.8[1:71])
spsecd
# show occurences with c higher sensitivity that d:
spsecd[spsecd$sens_c>=spsecd$sens_d,]
# small differences except for small n

# simple plots, specificities and sensitivities:
par(mar= c(bottom=4, left=3, top=3, right=.5)+.1)
par(mfrow=c(1,2))
plot(x=10:80, y=cr_bounds$pa_0.0[1:71], type="l", col="red", axes=FALSE,
      xlab="N", ylab="", main="Probability of true negative")
box()
axis(1, at=c(10,20,30,40,50,60,70,80))
axis(2, las=1)
points(x=10:80, y=cr_bounds$pb_0.0[1:71], col="green3", type="l")
points(x=10:80, y=cr_bounds$pc_0.0[1:71], col="blue", type="l")
points(x=10:80, y=cr_bounds$pd_0.0[1:71], col="orange", type="l")
# simple plot, sensitivities:
plot(x=10:80, y=1-cr_bounds$pa_0.8[1:71], type="l", col="red",
     axes=FALSE, xlab="N", ylab="", main="Sensitivity, shift 0.8")
box()
axis(1, at=c(10,20,30,40,50,60,70,80))
axis(2, las=1)
points(x=10:80, y=1-cr_bounds$pb_0.8[1:71], col="green3", type="l")
points(x=10:80, y=1-cr_bounds$pc_0.8[1:71], col="blue", type="l")
points(x=10:80, y=1-cr_bounds$pd_0.8[1:71], col="orange", type="l")
par(mfrow=c(1,1))
par(mar= c(bottom=5, left=4, top=4, right=2)+.1)

