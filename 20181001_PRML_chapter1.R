df = data.frame(theta= seq(-4,4,length.out = 301))
df$x  = dnorm(seq(-4,4,length.out = 301))
df$xl = dnorm(seq(-4,4,length.out = 301)*2 + 1)
df$x2 = dnorm((seq(-4,4,length.out = 301))^2)
df$x3 = dnorm((seq(-4,4,length.out = 301))^3)

df %>% tidyr::gather(key=g_y, value=dnorm,-theta) %>% 
  ggplot(aes(x=theta, y=dnorm, group=g_y))+
  geom_line()+
  facet_grid(g_y~.)




dice_l <- function(w,t=3,s=10,a=2,b=5){
  
  w # サイコロで1がでる確率
  v <- 1-w # サイコロで1以外がでる確率
  
  num <- prod(seq.int(s,s-t+1))
  den <- prod(seq.int(t:1))
  const <- num/den # 尤度の定数項，なくてもいい
  
  # 事後確率を求めるためにはサイコロすべての出目が分かる必要があるので，省略。
  res <- w^t * v^(s-t)* const * dbeta(w,a,b)# posterior ∝ likelihood * prior
  return(res)
}

dice_l(1/6)
dice_l(1/2)

# prior
ggplot(data=data.frame(x=c(0:1)),aes(x=x))+
  stat_function(fun=dbeta,args = list(shape1=2,shape2=5))

# posterior
ggplot(data =data.frame(x=c(0:1)),aes(x=x))+
  stat_function(fun=dice_l, args = list(t=4,s=10))



# 3D graphics

install.packages("mvtnorm")
library(mvtnorm)

dmvnorm(c(1,0),mean=rep(0,2), sigma = matrix(c(1,0.6,0.6,1),ncol=2))


mvnorm_df <- data.frame(x=seq(-4,4,length.out = 100) %>% rep(100),
                        y=seq(-4,4,length.out = 100) %>% rep(100) %>% sort())

mvnorm_df$z = dmvnorm(x=mvnorm_df,mean=rep(0,2), sigma = matrix(c(1,0.6,0.6,1),ncol=2))

mvnorm_df %>% ggplot(aes(x=x,y=y,z=z,fill=z))+
  geom_contour(aes(colour = stat(level)) ,binwidth=0.01)

# これを応用すれば，項目パラメタの等高線もかけるはず。
