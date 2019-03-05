library(tidyverse)
library(crossrun)
library(Rmpfr)

refresh <- F
n       <- 100
load('data/bounds.RData')

if(refresh) {
  load('data/crs100.RData')
  x <- list(
    x0 = map(crs100_0.0, asNumeric),
    x1 = map(crs100_1.0, asNumeric),
    x2 = map(crs100_2.0, asNumeric))
  rm(list = ls(pattern = 'crs100_'))
  write_rds(x, 'data/x.rds')
} else {
  x <- read_rds('data/x.rds')
}

crplot <- function(n = 12, shift = 0, labels = T, prop = F, round = 1) {
  ca <- bounds$ca[bounds$n == n]
  la <- bounds$la[bounds$n == n]
  pa <- bounds[n, paste0('pa_', shift, '.0')]
  cb <- bounds$cb[bounds$n == n]
  lb <- bounds$lb[bounds$n == n]
  pb <- bounds[n, paste0('pb_', shift, '.0')]
  cbord <- bounds$cbord[bounds$n == n]
  lbord <- bounds$lbord[bounds$n == n]
  pc <- bounds[n, paste0('pc_', shift, '.0')]
  
  m <- x[[paste0('x', shift)]][[paste0('pt', n)]]
  
  d <- m %>% 
    as_tibble(rownames = NA) %>% 
    rownames_to_column('C') %>% 
    gather('L', 'times', -C) %>% 
    mutate(L = as.numeric(L),
           C = as.numeric(C),
           p = times / sum(times),
           y = ifelse(prop * times, p, times),
           col = (times / max(times)) < 0.5 & (times / max(times)) > 0)
  
  p <- ggplot(d, aes(L, C, 
                     fill = times > 0,
                     alpha = times / max(times))) +
    geom_rect(aes(xmin = 0.55,
                  xmax = la - 0.55,
                  ymin = ca - 0.45,
                  ymax = max(C) + 0.45),
              colour = '#F8766D',
              fill = NA) +
    geom_rect(aes(xmin = 0.45,
                  xmax = lb - 0.45,
                  ymin = cb - 0.55,
                  ymax = max(C) + 0.55),
              colour = '#00BA38',
              fill = NA) +
    geom_rect(aes(xmin = lbord - 0.5,
                  xmax = lb - 0.5,
                  ymin = cb - 0.5,
                  ymax = cbord - 0.5),
              colour = '#619CFF',
              linetype = 5,
              fill = NA,
              na.rm = T) +
    geom_raster()
  
  if(labels) {
    p <- p +
      geom_text(aes(label = round(y, round),
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
          legend.position = 'none')
  # +
  #   labs(caption = paste('N =',
  #                        n,
  #                        '\nShift =',
  #                        shift,
  #                        '\nProportion captured by Anhøj rules =',
  #                        round(pa, 3)))
  plot(p)
}

crplot(10, labels = F)
crplot(11, labels = F)
crplot(15, labels = F)
crplot(16, labels = F)
crplot(19, labels = F)
crplot(47, labels = F)
