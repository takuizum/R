df = data.frame(theta= seq(-4,4,length.out = 301))
df$x  = dnorm(seq(-4,4,length.out = 301))
df$xl = dnorm(seq(-4,4,length.out = 301)*2 + 1)
df$x2 = dnorm((seq(-4,4,length.out = 301))^2)
df$x3 = dnorm((seq(-4,4,length.out = 301))^3)

df %>% tidyr::gather(key=g_y, value=dnorm,-theta) %>% 
  ggplot(aes(x=theta, y=dnorm, group=g_y))+
  geom_line()+
  facet_grid(g_y~.)
  
