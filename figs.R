# Setup ----
library(tidyverse)
library(crossrun)
load('data/bounds.RData')

refresh <- F

## Recalculate or get small list object for plotting LC box figures
if(refresh) {
  library(Rmpfr)
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

## Tall bounds data
limits <- bounds %>% 
  select(n:lb) %>% 
  gather('key', 'val', -n) %>% 
  separate(key, c('test', 'rule'), 1) %>% 
  mutate(test = fct_recode(test, `number of crossings` = 'c',
                           `longest run` = 'l'),
         rule = fct_recode(rule, anhoej = 'a',
                           `best box` = 'b'))
## Tall parameter value data
vals <- bounds %>% 
  select(-(ca:lbord)) %>%
  gather('key', 'val', -n) %>% 
  separate(key, c('test', 'shift'), '_') %>% 
  mutate(rule = substring(test, nchar(test)),
         test = substring(test, 1, nchar(test) - 1),
         shift = as.numeric(shift)) %>% 
  mutate(rule = fct_recode(rule, anhoej = 'a',
                           `best box` = 'b',
                           `cut box` = 'c'))

# Compare limits from anhoej and best box rules ----
ggplot(limits, aes(n, val, colour = rule)) +
  geom_line() +
  facet_wrap(~ test, scale = 'free_y') +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal()

# Plot power function ----
ggplot(filter(vals, test == 'p'), aes(n, 1 - val, colour = rule)) +
  geom_line() +
  facet_wrap(~ shift) +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Power function',
       y = 'Probability of signal',
       x = 'N')

## Specificity
ggplot(filter(vals, test == 'p', shift == 0), aes(n, val, colour = rule)) +
  geom_line() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Specificity',
       y = 'Probability of true negative',
       x = 'N')

# Plot likelihood ratios ----
## LR+
ggplot(filter(vals, test == 'lrpos'), aes(n, val, colour = rule)) +
  geom_line() +
  geom_hline(yintercept = 10) +
  facet_wrap(~ shift, ncol = 5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Positive likelihood ratios',
       y = 'LR+',
       x = 'N')

## LR-
ggplot(filter(vals, test == 'lrneg'), aes(n, val, colour = rule)) +
  geom_line() +
  geom_hline(yintercept = 0.1) +
  facet_wrap(~ shift, ncol = 5) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(20, 100, by = 20)) +
  theme_minimal() +
  labs(title = 'Negative likelihood ratios',
       y = 'LR-',
       x = 'N')

# Function to plot LC box figures ----
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
    geom_raster() +
    geom_rect(aes(xmin = 0.5,                       # Anhoej box
                  xmax = la + 0.5,
                  ymin = ca - 0.5,
                  ymax = max(C) + 0.5),
              colour = '#F8766D',
              fill = NA) +
    geom_rect(aes(xmin = 0.5,                       # Best box
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = max(C) + 0.5),
              linetype = 2,
              colour = '#00BA38',
              fill = NA) +
    geom_rect(aes(xmin = lbord + 0.5,               # Cut box, horizontal
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = cb + 0.5),
              colour = '#619CFF',
              fill = NA,
              linetype = 4,
              na.rm = T) +
    geom_rect(aes(xmin = lb - 0.5,                  # Cut box, vertical
                  xmax = lb + 0.5,
                  ymin = cb - 0.5,
                  ymax = cbord - 0.5),
              colour = '#619CFF', #619CFF',
              linetype = 4,
              fill = NA,
              na.rm = T)
  
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

# crplot(10, labels = F)
# crplot(11, labels = F)
# crplot(15, labels = F)
crplot(16, labels = F)
# crplot(19, labels = F)
# crplot(24, labels = F)
# crplot(47, labels = F)
# crplot(52, labels = F)

