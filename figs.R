library(tidyverse)
load('crossrunbox1.Rdata')

x <- bounds %>% 
  transmute(n,
         c_anhoej = ca,
         l_anhoej = la,
         c_bestbox = cb,
         l_bestbox = lb,
         c_cutbox = ifelse(is.na(cbord), c_bestbox, cbord),
         l_cutbox = ifelse(is.na(lbord), l_bestbox, lbord),
         specificity_anhoej = pa.0,
         specificity_bestbox = pb.0,
         specificity_cutbox = pc.0,
         lrposa.1,
         lrposb.1,
         lrposc.1,
         lrposa.2,
         lrposb.2,
         lrposc.2,
         lrnega.1,
         lrnegb.1,
         lrnegc.1,
         lrnega.2,
         lrnegb.2,
         lrnegc.2)
         
ggplot(x[1:30,], aes(n, lrposa.2)) +
  geom_line() +
  geom_line(aes(y = lrnega.2)) +
  scale_y_log10(breaks = c(0.01, 0.1, 0, 10)) +
  scale_x_continuous(breaks = seq(10, 100, by = 10))

ggplot(x[1:30,], aes(n, lrposb.2)) +
  geom_line() +
  geom_line(aes(y = lrnegb.2)) +
  scale_y_log10(breaks = c(0.01, 0.1, 0, 10))

ggplot(x[1:30,], aes(n, lrposc.2)) +
  geom_line() +
  geom_line(aes(y = lrnegc.2)) +
  scale_y_log10(breaks = c(0.01, 0.1, 0, 10))

s <- x %>% 
  select(n, specificity_anhoej:specificity_cutbox) %>% 
  gather('key', 'val', -n)

ggplot(s, aes(n, val)) +
  geom_line() +
  facet_wrap(~ key)

lrp <- x %>% 
  select(n, lrposa.2:lrposc.2) %>% 
  gather('key', 'val', -n)

ggplot(lrp, aes(n, val)) +
  geom_line() +
  facet_wrap(~ key) +
  scale_y_log10(breaks = c(1, 10, 100))

lrn <- x %>% 
  select(n, lrnega.2:lrnegc.2) %>% 
  gather('key', 'val', -n)

ggplot(lrn, aes(n, val)) +
  geom_line() +
  facet_wrap(~ key) +
  scale_y_log10(breaks = c(0, 0.1, 0.01))
