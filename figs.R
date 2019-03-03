library(tidyverse)
bounds <- readRDS('bounds.rds')

# Compare limits from anhoej and best box rules
limits <- bounds %>% 
  select(n:lb) %>% 
  gather('key', 'val', -n) %>% 
  separate(key, c('test', 'rule'), 1) %>% 
  mutate(test = fct_recode(test, `number of crossings` = 'c',
                           `longest run` = 'l'),
         rule = fct_recode(rule, anhoej = 'a',
                           `best box` = 'b'))

ggplot(limits, aes(n, val, colour = rule)) +
  geom_line() +
  facet_wrap(~ test, scale = 'free_y')

# Wide to tall data
tall <- bounds %>% 
  select(-(ca:lbord)) %>%
  gather('key', 'val', -n) %>% 
  separate(key, c('test', 'shift'), '_') %>% 
  mutate(rule = substring(test, nchar(test)),
         test = substring(test, 1, nchar(test) - 1),
         shift = as.numeric(shift)) %>% 
  mutate(rule = fct_recode(rule, anhoej = 'a',
                           `best box` = 'b',
                           `cut box` = 'c'))


# Plot power function
ggplot(filter(tall, test == 'p'), aes(n, 1 - val, colour = rule)) +
  geom_line() +
  facet_wrap(~ shift) +
  theme_minimal() +
  labs(title = 'Power function',
       y = 'Probability of signal',
       x = 'N')

# Plot LR+
ggplot(filter(tall, test == 'lrpos'), aes(n, val, colour = rule)) +
  geom_line() +
  geom_hline(yintercept = 10) +
  facet_wrap(~ shift) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = 'Positive likelihood ratio',
       y = 'LR+',
       x = 'N')

# Plot LR-
ggplot(filter(tall, test == 'lrneg'), aes(n, val, colour = rule)) +
  geom_line() +
  geom_hline(yintercept = 0.1) +
  facet_wrap(~ shift) +
  scale_y_log10() +
  theme_minimal() +
  labs(title = 'Negative likelihood ratio',
       y = 'LR-',
       x = 'N')
