ggplot(data=data.frame((x=c(0:1))), aes(x=x))+
  stat_function(fun=dbinom, args=list(x=7,size=10))

set.seed(123)
rbinom(10000, size = 30, prob = 0.8) %>% data.frame() %>%
  ggplot(aes(x=.))+
  geom_histogram(bins = 16)
set.seed(123)
rbinom(10000, size = 30, prob = 0.8) %>% scale() %>%  data.frame() %>% 
  ggplot(aes(x=.))+
  geom_histogram(bins = 15)



rmultinom(10,size=c(1:6),prob=rep(1/6,6))



# normal ogive model
#install.packages("latex2exp")
library(latex2exp)
library(ggrepel)
library(tidyverse)
library(devEMF)
library(irtfun2)
library(directlabels)


#1-1
emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-1_normal1.emf", width = 7, height = 3)
plot(dnorm, -4, 4, ann=F, axes=F)
n  <- 100
xs <- seq(-1.5, 4, length=n)
ys <- dnorm(xs)
polygon(c(xs,rev(xs)), c(rep(0,n),rev(ys)), col="gray")
abline(h=0)
abline(v=0, col="blue")
dev.off()
#1-2
emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-1_normal2.emf", width = 7, height = 3)
plot(dnorm, -4, 4, ann=F, axes=F)
n  <- 100
xs <- seq(0.2, 4, length=n)
ys <- dnorm(xs)
polygon(c(xs,rev(xs)), c(rep(0,n),rev(ys)), col="gray")
abline(h=0)
abline(v=0, col="blue")
dev.off()
#1-3
emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-1_normal3.emf", width = 7, height = 3)
plot(dnorm, -4, 4, ann=F, axes=F)
n  <- 100
xs <- seq(2, 4, length=n)
ys <- dnorm(xs)
polygon(c(xs,rev(xs)), c(rep(0,n),rev(ys)), col="gray")
abline(h=0)
abline(v=0, col="blue")
dev.off()

# 項目反応モデルの理論的根拠の図
icc <- function(theta,a,b,c=0,D){
  e <- exp(-D*a*(theta-b))
  p <- c+(1-c)/(1+e)
}

nom_g <- 
  ggplot(data=data.frame((x=c(-4:4))), aes(x=x))+
  stat_function(fun=pnorm, args=list(mean=0,sd=1), aes(colour="Normal_Ogive"))+
  stat_function(fun=icc, args = list(a=1,b=0,D=1.702), aes(colour="D=1.702"))+
  stat_function(fun=icc, args = list(a=1,b=0,D=1.749), aes(colour="D=1.749"))+
  labs(colour="type")+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-2_NOM.emf", width = 7, height = 5)
nom_g
dev.off()


# ICC

icc_1pl_g <- ggplot(data = data.frame(x=c(-4:4)), aes(x=x))+
  stat_function(fun=icc, args = list(a=1,b=0,D=1.702), aes(colour="b=0"))+
  stat_function(fun=icc, args = list(a=1,b=-2,D=1.702), aes(colour="b=-2"))+
  stat_function(fun=icc, args = list(a=1,b=1.5,D=1.702), aes(colour="b=1.5"))+
  labs(colour="difficulty")+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-3_icc1pl.emf", width = 7, height = 5)
icc_1pl_g
dev.off()
  
icc_2pl_g <- ggplot(data = data.frame(x=c(-4:4)), aes(x=x))+
  stat_function(fun=icc, args = list(a=0.5,b=0,D=1.702), aes(colour="a=0.5, b=0"))+
  stat_function(fun=icc, args = list(a=1.2,b=-2,D=1.702), aes(colour="a=1.2, b=-2"))+
  stat_function(fun=icc, args = list(a=2,b=1.5,D=1.702), aes(colour="a=2, b=1.5"))+
  labs(colour="discrimination\n & difficulty")+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-4_icc2pl.emf", width = 7, height = 5)
icc_2pl_g
dev.off()

icc_3pl_g <- ggplot(data = data.frame(x=c(-4:4)), aes(x=x))+
  stat_function(fun=icc, args = list(a=0.5,b=0,c=0.2,D=1.702), aes(colour="a=0.5, b=0, c=0.2"))+
  stat_function(fun=icc, args = list(a=1.2,b=-2,c=0.1,D=1.702), aes(colour="a=1.2, b=-2, c=0.1"))+
  stat_function(fun=icc, args = list(a=2,b=1.5,c=0.1,D=1.702), aes(colour="a=2, b=1.5, c=0.1"))+
  ylim(0,1)+
  labs(colour="discrimination\n & difficulty\n & guessing")+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-5_icc3pl.emf", width = 7, height = 5)
icc_3pl_g
dev.off()


# GPCM

gpcm_sub <- function(theta, a, b, D, k){
  K <- length(b)
  G <- rep(1,K+1)
  for(v in 1:K) G[v+1] <- exp(sum(D*a*(theta-b[1:v])))
  p <- G[k+1]/sum(G)
  p
}

gpcm <- function(theta, a, b, D, k){
  apply(as.matrix(theta), 1, gpcm_sub, a=a,b=b,D=D,k=k)
}

tp <- c(-1.5, 0, 1)

gpcm_g <- ggplot(data = data.frame(x=c(-4:4)), aes(x=x))+
  stat_function(fun=gpcm, args = list(a=1.5, b=tp, D=1.702, k=0), aes(colour="category0"))+
  stat_function(fun=gpcm, args = list(a=1.5, b=tp, D=1.702, k=1), aes(colour="category1"))+
  stat_function(fun=gpcm, args = list(a=1.5, b=tp, D=1.702, k=2), aes(colour="category2"))+
  stat_function(fun=gpcm, args = list(a=1.5, b=tp, D=1.702, k=3), aes(colour="category3"))+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"), colour="Category")

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-6_gpcm.emf", width = 7, height = 5)
gpcm_g
dev.off()


# GIRT model
gptheta <- function(theta,phi,a,b,D){
  A <- sqrt(1+phi^2*a^2)
  e <- exp(-D*a/A*(theta-b))
  p <- 1/(1+e)
  p
}

girt_g <- ggplot(data = data.frame(x=c(-4:4)), aes(x=x))+
  stat_function(fun=gptheta, args = list(a=1,b=0,D=1.702, phi=0.1), aes(colour="phi=0.1"))+
  stat_function(fun=gptheta, args = list(a=1,b=0,D=1.702, phi=1), aes(colour="phi=1"))+
  stat_function(fun=gptheta, args = list(a=1,b=0,D=1.702, phi=2), aes(colour="phi=2"))+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"), colour=TeX("$\\phi$"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/1-7_girt-icc.emf", width = 7, height = 5)
girt_g
dev.off()


# Newton-Raphson & Fisher's Scoring


# P(theta) in two-parameter logisticmodel
ptheta <- function(theta,a,b,c,D=1.702){
  c+(1-c)/(1+exp(-D*a*(theta-b)))
}

iif <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  I <- D^2*a^2*(1-p)*(p-c)^2/((1-c)^2*p)
}

# 対数尤度
LL <- function(u,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  sum(u*log(p)+(1-u)*log(1-p),na.rm = T)
}

# 一階偏微分
fpd <- function(xi,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  D*sum(a*(xi - p)*(p-c)/(p*(1-c)), na.rm = T)
}

# 二階偏微分
spd <- function(xi,theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T)
}

# テスト情報量（尤度関数の二階偏微分の負の期待値）
pitheta <- function(theta,a,b,c,D){
  p <- ptheta(theta,a,b,c,D)
  I <- D^2*sum(a^2*(1-p)*(p-c)^2/((1-c)^2*p), na.rm = T)
}

line_ic <- function(y,a,x,X){
  #傾きとx,yの値から切片を求めて直線の値yを再計算する関数
  b <- y-(a*x)
  a*X+b
}

line_ic2 <- function(y,a,x,X){
  #傾きとx,yの値から切片を求めて直線の値yを再計算する関数
  a <- -a
  b <- y-(a*x)
  a*X+b
} 

nr_fun <- function(t0,fp,fpp){
  # ニュートンラフソンの更新した値を計算する関数
  t0-fp/fpp
}

fs_fun <- function(t0,fp,I){
  # フィッシャースコアリングの更新した値を計算する関数
  t0+fp/I
}


# response pattern"
set.seed(123)
u <- sample(c(0,1),30, replace = T)
# item parameters
# 2PLM
a <- rlnorm(30,meanlog = -0.5,sdlog = 0.5)
b <- rnorm(30)
c <- rep(0,30)
#c <- rbeta(30,2,10)

dat_L <- data.frame(x=seq(-4,4,length.out = 101))
dat_L$LL <- apply(matrix(dat_L$x), 1, LL, u=u,a=a,b=b,c=c,D=1.702) 
dat_L$fpd <- apply(matrix(dat_L$x), 1, fpd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$spd <- apply(matrix(dat_L$x), 1, spd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$I <- -apply(matrix(dat_L$x), 1, pitheta, a=a,b=b,c=c,D=1.702) 

LL_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=LL))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL(\\theta)$"))

fpd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=fpd))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"))

spd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=spd))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL''(\\theta)$"))

info_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=I))+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-1_loglikelihoodF.emf", width = 7, height = 5)
LL_g
dev.off()

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-2_fpd.emf", width = 7, height = 5)
fpd_g
dev.off()

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-3_spd.emf", width = 7, height = 5)
spd_g
dev.off()

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-4_minusinfo.emf", width = 7, height = 5)
info_g
dev.off()


# まずは-2から
NR_g <- fpd_g+
  stat_function(fun=line_ic, args = list(y=fpd(u,-2,a,b,c,D=1.702), a=spd(u,-2,a,b,c,D=1.702), x=-2), aes(colour=1))
t2 <- nr_fun(-2,fpd(u,-2,a,b,c,D=1.702),spd(u,-2,a,b,c,D=1.702))

# 自動作図
NR_g <- fpd_g
t0 <- -2
t <- 0
conv <- T
while(conv){
  t <- t+1
  slope <- spd(u,t0,a,b,c,D=1.702)
  NR_g <- NR_g+
    stat_function(fun=line_ic, args = list(y=fpd(u,t0,a,b,c,D=1.702), a=slope, x=t0),
                  aes_q(colour=sprintf("iterations=%d, theta=%.3f",t,t0)))
  t1 <- nr_fun(t0,fpd(u,t0,a,b,c,D=1.702),spd(u,t0,a,b,c,D=1.702))
  if(abs(t1-t0)<0.001) conv <- FALSE
  t0 <- t1
}
NR_g <- NR_g+labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"), colour=TeX("iterations"))

FS_g <- fpd_g
t0 <- -2
t <- 0
conv <- T
while(conv){
  t <- t+1
  I <- pitheta(t0,a,b,c,D=1.702)
  FS_g <- FS_g+
    stat_function(fun=line_ic2, args = list(y=fpd(u,t0,a,b,c,D=1.702), a=I, x=t0),
                  aes_q(colour=sprintf("iterations=%d, theta=%.3f",t,t0)))
  t1 <- fs_fun(t0,fpd(u,t0,a,b,c,D=1.702),pitheta(t0,a,b,c,D=1.702))
  if(abs(t1-t0)<0.001) conv <- FALSE
  t0 <- t1
}
FS_g <- FS_g+labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"), colour=TeX("iterations"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-5_NR.emf", width = 7, height = 5)
NR_g
dev.off()

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-6_FS.emf", width = 7, height = 5)
FS_g
dev.off()


# 3PLM
set.seed(123)
a <- rlnorm(30,meanlog = -0.5,sdlog = 0.5)
b <- rnorm(30)
c <- rbeta(30,2,10)

dat_L <- data.frame(x=seq(-4,4,length.out = 101))
dat_L$LL <- apply(matrix(dat_L$x), 1, LL, u=u,a=a,b=b,c=c,D=1.702) 
dat_L$fpd <- apply(matrix(dat_L$x), 1, fpd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$spd <- apply(matrix(dat_L$x), 1, spd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$I <- -apply(matrix(dat_L$x), 1, pitheta, a=a,b=b,c=c,D=1.702) 

LL_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=LL))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL(\\theta)$"))

fpd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=fpd))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"))

spd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=spd))+
  labs(x=TeX("$\\theta$"), y=TeX("$lnL''(\\theta)$"))

info_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=I))+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-7_loglikelihoodF_3PL.emf", width = 7, height = 5)
LL_g
dev.off()

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-8_fpd_3PL.emf", width = 7, height = 5)
fpd_g
dev.off()

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-9_spd_3PL.emf", width = 7, height = 5)
spd_g
dev.off()

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-10_minusinfo_3PL.emf", width = 7, height = 5)
info_g
dev.off()

# 自動作図
NR_g <- fpd_g
t0 <- 0
t <- 0
conv <- T
while(conv){
  t <- t+1
  slope <- spd(u,t0,a,b,c,D=1.702)
  NR_g <- NR_g+
    stat_function(fun=line_ic, args = list(y=fpd(u,t0,a,b,c,D=1.702), a=slope, x=t0),
                  aes_q(colour=sprintf("iterations=%d, theta=%.3f",t,t0)))
  t1 <- nr_fun(t0,fpd(u,t0,a,b,c,D=1.702),spd(u,t0,a,b,c,D=1.702))
  if(abs(t1-t0)<0.001) conv <- FALSE
  t0 <- t1
}
NR_g <- NR_g+labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"), colour=TeX("iterations"))

FS_g <- fpd_g
t0 <- 0
t <- 0
conv <- T
while(conv){
  t <- t+1
  I <- pitheta(t0,a,b,c,D=1.702)
  FS_g <- FS_g+
    stat_function(fun=line_ic2, args = list(y=fpd(u,t0,a,b,c,D=1.702), a=I, x=t0),
                  aes_q(colour=sprintf("iterations=%d, theta=%.3f",t,t0)))
  t1 <- fs_fun(t0,fpd(u,t0,a,b,c,D=1.702),pitheta(t0,a,b,c,D=1.702))
  if(abs(t1-t0)<0.001) conv <- FALSE
  t0 <- t1
}
FS_g <- FS_g+labs(x=TeX("$\\theta$"), y=TeX("$lnL'(\\theta)$"), colour=TeX("iterations"))

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-11_NR_3PL.emf", width = 7, height = 5)
NR_g
dev.off()

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-12_FS_3PL.emf", width = 7, height = 5)
FS_g
dev.off()


# MAP

# 尤度関数
LL_b <- function(u,theta,a,b,c,mu,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  sum(log(p)*u+log(1-p)*(1-u)-0.5*((theta-mu)/sigma)^2,na.rm = T)
}

# 一階偏微分
fpdLPD <- function(xi,theta,a,b,c,mu,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  D*sum(a*(xi - p)*(p-c)/(p*(1-c)), na.rm = T)-1/sigma^2*(theta-mu)
}

# 二階偏微分
spdLPD <- function(xi,theta,a,b,c,sigma,D){
  p <- ptheta(theta,a,b,c,D)
  D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T) - 1/sigma
}

set.seed(123)
a <- rlnorm(30,meanlog = -0.5, sdlog = 0.5)
b <- rnorm(30)
c <- rbeta(30,2,10)

dat_L <- data.frame(x=seq(-4,4,length.out = 101))
dat_L$LL <- apply(matrix(dat_L$x), 1, LL, u=u,a=a,b=b,c=c,D=1.702) 
dat_L$fpd <- apply(matrix(dat_L$x), 1, fpd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$spd <- apply(matrix(dat_L$x), 1, spd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$LL_B <- apply(matrix(dat_L$x), 1, LL_b, u=u,a=a,b=b,c=c,D=1.702,mu=0,sigma=1) 
dat_L$fpdLPD <- apply(matrix(dat_L$x), 1, fpdLPD, xi=u,a=a,b=b,c=c,D=1.702,mu=0,sigma=1) 
dat_L$spdLPD <- apply(matrix(dat_L$x), 1, spdLPD, xi=u,a=a,b=b,c=c,D=1.702,sigma=1) 

MAP_LL_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=LL,colour="MLE"))+
  geom_line(aes(y=LL_B, colour="MAP"))+
  labs(x=TeX("$\\theta$"),y="",colour="mathod")

MAP_fpd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=fpd,colour="MLE"))+
  geom_line(aes(y=fpdLPD, colour="MAP"))+
  labs(x=TeX("$\\theta$"),y="",colour="mathod")

MAP_spd_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=spd,colour="MLE"))+
  geom_line(aes(y=spdLPD, colour="MAP"))+
  labs(x=TeX("$\\theta$"),y="",colour="mathod")

# 事前分布のパラメタを変化させた場合の尤度方程式の変化
dat_L <- data.frame(x=seq(-4,4,length.out = 101))
dat_L$fpd <- apply(matrix(dat_L$x), 1, fpd, xi=u,a=a,b=b,c=c,D=1.702) 
dat_L$fpdLPD1 <- apply(matrix(dat_L$x), 1, fpdLPD, xi=u,a=a,b=b,c=c,D=1.702,mu=0,sigma=1) 
dat_L$fpdLPD2 <- apply(matrix(dat_L$x), 1, fpdLPD, xi=u,a=a,b=b,c=c,D=1.702,mu=5,sigma=1) 
dat_L$fpdLPD3 <- apply(matrix(dat_L$x), 1, fpdLPD, xi=u,a=a,b=b,c=c,D=1.702,mu=0,sigma=0.8) 

MAP_fpd2_g <- ggplot(dat_L, aes(x=x))+
  geom_line(aes(y=fpd,colour="MLE"))+
  geom_line(aes(y=fpdLPD1, colour="N(0,1)"))+
  geom_line(aes(y=fpdLPD2, colour="N(5,1)"))+
  geom_line(aes(y=fpdLPD3, colour="N(0,0.64)"))+
  labs(x=TeX("$\\theta$"),y=TeX("$ln\\L'(\\theta)+lnp'(\\theta)$"),colour="condition")

emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/2-13_MAP_fpd2.emf", width = 7, height = 5)
MAP_fpd2_g
dev.off()


# 実験結果のプロット(Vertical Scalingプロジェクト上で実行のこと)




#-------------------------------------------------------#
#edit "japanese test" data
#-------------------------------------------------------#

e4 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-e4-j.csv")
ID　<- e4$個人番号
pop.num <- rep(1,nrow(e4))
e4_j <- cbind(ID,pop.num,e4[,-1])
e5 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-e5-j.csv")
ID　<- e5$個人番号
pop.num <- rep(2,nrow(e5))
e5_j <- cbind(ID,pop.num,e5[,-1])
e6 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-e6-j.csv")
ID　<- e6$個人番号
pop.num <- rep(3,nrow(e6))
e6_j <- cbind(ID,pop.num,e6[,-1])
j1 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-j1-j.csv")
ID　<- j1$個人番号
pop.num <- rep(4,nrow(j1))
j1_j <- cbind(ID,pop.num,j1[,-1])
j2 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-j2-j.csv")
ID　<- j2$個人番号
pop.num <- rep(5,nrow(j2))
j2_j <- cbind(ID,pop.num,j2[,-1])


full_data_j <- e4_j %>% 
  dplyr::full_join(e5_j) %>% 
  dplyr::full_join(e6_j) %>% 
  dplyr::full_join(j1_j) %>% 
  dplyr::full_join(j2_j)



#-------------------------------------------------------#
#edit math data
#-------------------------------------------------------#

e4 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-e4-m.csv")
ID　<- e4$個人番号
pop.num <- rep(1,nrow(e4))
e4_m <- cbind(ID,pop.num,e4[,-1])
e5 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-e5-m.csv")
ID　<- e5$個人番号
pop.num <- rep(2,nrow(e5))
e5_m <- cbind(ID,pop.num,e5[,-1])
e6 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-e6-m.csv")
ID　<- e6$個人番号
pop.num <- rep(3,nrow(e6))
e6_m <- cbind(ID,pop.num,e6[,-1])
j1 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-j1-m.csv")
ID　<- j1$個人番号
pop.num <- rep(4,nrow(j1))
j1_m <- cbind(ID,pop.num,j1[,-1])
j2 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/00_edited data/data/h29-j2-m.csv")
ID　<- j2$個人番号
pop.num <- rep(5,nrow(j2))
j2_m <- cbind(ID,pop.num,j2[,-1])


full_data_m <- e4_m %>% 
  dplyr::full_join(e5_m) %>% 
  dplyr::full_join(e6_m) %>% 
  dplyr::full_join(j1_m) %>% 
  dplyr::full_join(j2_m)


#check test unidimensionality
library(psych)

#japanese

eigen_j <- data.frame(component_number=c(1:30), 
                      G1=eigen(tetrachoric(e4_j[,c(-1,-2)])$rho)$values,
                      G2=eigen(tetrachoric(e5_j[,c(-1,-2)])$rho)$values,
                      G3=eigen(tetrachoric(e6_j[,c(-1,-2)])$rho)$values,
                      G4=eigen(tetrachoric(j1_j[,c(-1,-2)])$rho)$values,
                      G5=eigen(tetrachoric(j2_j[,c(-1,-2)])$rho)$values)

g_eigen_j <- eigen_j %>% tidyr::gather(key=grade, value=Eigen_value, -component_number) %>% 
  ggplot(aes(x=component_number, y=Eigen_value, colour=grade))+
  geom_point()+
  geom_line()

eigen_j_st <- data.frame(component_number=c(1:10), 
                         Eigen_value=eigen(tetrachoric(
                           full_data_j %>% dplyr::select(starts_with("st"))
                         )$rho)$values)

g_eigen_j_st <- eigen_j_st %>% 
  ggplot(aes(x=component_number, y=Eigen_value))+
  geom_point()+
  geom_line()

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-1_eigen_j.emf", width = 7, height = 5)
g_eigen_j
dev.off()


devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-2_eigen_j_st.emf", width = 7, height = 5)
g_eigen_j_st
dev.off()

#math
eigen_m <- data.frame(component_number=c(1:30), 
                      G1=eigen(tetrachoric(e4_m[,c(-1,-2)])$rho)$values,
                      G2=c(eigen(tetrachoric(e5_m[,c(-1,-2)])$rho)$values,0),
                      G3=eigen(tetrachoric(e6_m[,c(-1,-2)])$rho)$values,
                      G4=eigen(tetrachoric(j1_m[,c(-1,-2)])$rho)$values,
                      G5=eigen(tetrachoric(j2_m[,c(-1,-2)])$rho)$values)

g_eigen_m <- eigen_m %>% tidyr::gather(key=grade, value=Eigen_value, -component_number) %>% 
  ggplot(aes(x=component_number, y=Eigen_value, colour=grade))+
  geom_point()+
  geom_line()

eigen_m_st <- data.frame(component_number=c(1:10), 
                         Eigen_value=eigen(tetrachoric(
                           full_data_m %>% dplyr::select(starts_with("st"))
                         )$rho)$values)

g_eigen_m_st <- eigen_m_st %>% 
  ggplot(aes(x=component_number, y=Eigen_value))+
  geom_point()+
  geom_line()


devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-3_eigen_m.emf", width = 7, height = 5)
g_eigen_m
dev.off()

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-4_eigen_m_st.emf", width = 7, height = 5)
g_eigen_m_st
dev.off()



# Item Analysis
res_con_j <- full_data_j %>% estip(fc=3, bg=3, ng=5, min_a = 0, maxabs_b = 20)
res_con_j$para
res_con_m <- full_data_m %>% estip(fc=3, bg=3, ng=5, min_a = 0, maxabs_b = 20)
res_con_m$para

res_eap_j <- full_data_j %>% estheta(param=res_con_j$para, sigma=3)
res_eap_m <- full_data_m %>% estheta(param=res_con_m$para, sigma=3)


Q3_j <- Q3_stat(x=full_data_j, para = res_con_j$para, theta=res_eap_j$res$EAP, abs_value = 0.2)
Q3_j$df[isTRUE(Q3_j$df$abs)]
View(Q3_j$matrix)
Q3_m <- Q3_stat(x=full_data_m, para = res_con_m$para, theta=res_eap_m$res$EAP, abs_value = 0.2)
Q3_m$df[isTRUE(Q3_m$df$abs)]
View(Q3_m$matrix)


# item fit

fit_j <- ifind3(x=full_data_j, para = res_con_j$para, theta=res_eap_j$res$EAP)
fit_m <- ifind3(x=full_data_m, para = res_con_m$para, theta=res_eap_m$res$EAP)


fit_j2 <- ifind(x=full_data_j, para = res_con_j$para, theta=res_eap_j$res$EAP, line_size = 0.5)
fit_j2$ggplot
fit_m2 <- ifind(x=full_data_m, para = res_con_m$para, theta=res_eap_m$res$EAP, line_size = 0.5)
fit_m2$ggplot
fit_m2$X2

write.csv(fit_j, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-9-1_itemfit_temp_j.csv", quote = F, row.names = F)
write.csv(fit_m, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-9-2_itemfit_temp_m.csv", quote = F, row.names = F)


ggsave(path = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/", filename = "4-21_fit_graph_j.png")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-21_fit_graph_j.emf", width = 5, height = 9)
fit_j2$ggplot+
  theme(axis.text = element_text(size=7),
        strip.text = element_text(size=7),
        strip.background = element_blank()) # facet先のタイトルの文字の大きさを変更する。
dev.off()

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-22_fit_graph_m.emf", width = 5, height = 9)
fit_m2$ggplot+
  theme(axis.text = element_text(size=7),
        strip.text = element_text(size=7),
        strip.background = element_blank())
dev.off()

write.csv(res_con_j$para, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-3_itempara_temp_j.csv", quote = F, row.names = F)
write.csv(res_con_m$para, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-3_itempara_temp_m.csv", quote = F, row.names = F)


# re estimate

# prior dist setting
# theta ~ Normal(0,1) # multigroup model

# japanese
res_j <- full_data_j %>% estip(fc=3, bg=3, ng=5, rm_list = c("st101","st108"), min_a = 0)

res_j$para
res_j$SE

write.csv(res_j$para, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-12_itempara_j.csv", quote = F, row.names = F)
write.csv(res_j$SE, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-13_itempara_se_j.csv", quote = F, row.names = F)

eap_j <- full_data_j %>% estheta(param=res_j$para, mintheta = -4, maxtheta = 4, sigma = 3)

fit_j3 <- ifind(full_data_j, res_j$para, eap_j$res$EAP)
fit_j3$ggplot
fit_j3$X2
fit_j4 <- ifind3(full_data_j, res_j$para, eap_j$res$EAP)
write.csv(fit_j4, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-9-3_itemfit_j.csv", quote = F, row.names = F)


Q3_j2 <- Q3_stat(full_data_j, res_j$para, eap_j$res$EAP, abs_value = 0.3)
Q3_j2$df[Q3_j2$df$abs %>% isTRUE(),]

j <- res_j$para$a != 0 # key
res_j$para$level <- res_j$para$Item %>% str_extract("[a-z]{1,2}")
#ICC
d_icc_j <- apply(matrix(seq(-6,6,length.out = 101)), 1, ptheta, 
                 a=res_j$para$a[j], b=res_j$para$b[j], c=res_j$para$c[j], D=1.702) %>% t() %>% data.frame()
colnames(d_icc_j) <- res_j$para$Item[j]
d_icc_j$theta <- seq(-6,6,length.out = 101)
g_icc_j <- d_icc_j %>% tidyr::gather(key=Item, value=prob, -theta) %>% 
  ggplot(aes(x=theta, y=prob, group=Item, colour=Item))+
  geom_line()+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"), colour=TeX("Item"))
#IIF
d_iif_j <- apply(matrix(seq(-6,6,length.out = 101)), 1, iif, 
                 a=res_j$para$a[j], b=res_j$para$b[j], c=res_j$para$c[j], D=1.702) %>% t() %>% data.frame()
colnames(d_iif_j) <- res_j$para$Item[j]
d_iif_j$theta <- seq(-6,6,length.out = 101)
g_iif_j <- d_iif_j %>% tidyr::gather(key=Item, value=prob, -theta) %>% 
  ggplot(aes(x=theta, y=prob, group=Item, colour=Item))+
  geom_line()+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"), colour=TeX("Item"))

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-23_icc_j.emf", width = 7, height = 5)
g_icc_j
dev.off()

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-24_iif_j.emf", width = 7, height = 5)
g_iif_j+ylim(0,3.8)
dev.off()

#scater
g_item_scat_j <- res_j$para[j,] %>% ggplot(aes(x=b, y=a, colour=Item))+
  geom_point()+xlim(-7,5)+ylim(0,2.5)+
  labs(x="困難度", y="識別力", colour="Item")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-27_item_scatter_j.emf", width = 7, height = 5)
g_item_scat_j
dev.off()

#scatter 2
g_item_scat_j2 <- res_j$para[j,] %>% ggplot(aes(x=b, y=a, colour=level))+
  geom_point()+xlim(-7,5)+ylim(0,2.5)+
  geom_text_repel(aes(label=Item))+
  labs(x="困難度", y="識別力", colour="Item")+
  facet_wrap(~level, ncol = 1)+
  theme(legend.position = 'none')

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-31_item_scat_lv_j.emf", width = 7, height = 9)
g_item_scat_j2
dev.off()

# TIF
d_tif_j <- apply(matrix(seq(-6,6,length.out = 101)), 1, tif, 
                 a=res_j$para$a[j], b=res_j$para$b[j], c=res_j$para$c[j], D=1.702) %>% data.frame()
d_tif_j$theta <- seq(-6,6,length.out = 101)
g_tif_j <- d_tif_j %>%
  ggplot(aes(x=theta, y=.))+
  geom_line()+ylim(0,21)+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))+
  theme(legend.position = 'none')

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-29_tif_j.emf", width = 7, height = 9)
g_tif_j
dev.off()

# TIF of each levels
res_j$para$level <- res_j$para$Item %>% str_extract("[a-z]{1,2}")
d_tif_lv_j <- data.frame(a=tif_auto(res_j$para[res_j$para$level=="a",], max = 10, min = -10),
                         b=tif_auto(res_j$para[res_j$para$level=="b",], max = 10, min = -10),
                         c=tif_auto(res_j$para[res_j$para$level=="c",], max = 10, min = -10),
                         d=tif_auto(res_j$para[res_j$para$level=="d",], max = 10, min = -10),
                         e=tif_auto(res_j$para[res_j$para$level=="e",], max = 10, min = -10),
                         f=tif_auto(res_j$para[res_j$para$level=="f",], max = 10, min = -10),
                         st=tif_auto(res_j$para[res_j$para$level=="st",], max = 10, min = -10),
                         theta=seq(-6,6,length.out = 301))
g_tif_lv_j <- d_tif_lv_j %>% tidyr::gather(key=levels, value=tif, -theta) %>% 
  ggplot(aes(x=theta, y=tif, colour=levels))+
  geom_line(aes(linetype=levels), size=1.5)+#ylim(0,6)+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))+
  annotate("text", x=c(-7,-4,-1,-2.5,-1,3,-2)*0.4, y=apply(d_tif_lv_j,2,max)[-8]+c(0,0.1,0,0.05,0,0,-0.2), 
           label=c("a","b","c","d","e","f","st"), colour=scales::hue_pal()(7), size=5, fontface="bold")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-32_tif_lv_j.emf", width = 7, height = 5)
g_tif_lv_j
dev.off()


# population distribution
d_pdist_j <- res_j$population_dist %>% data.frame() 
colnames(d_pdist_j) <- c("theta","G1","G2","G3","G4","G5")
g_pdist_j <- d_pdist_j %>% tidyr::gather(key=grade, value=Probability, -theta) %>% 
  ggplot(aes(x=theta, y=Probability, colour=grade))+
  geom_line()+
  annotate("text", x=c(-3,-0.5,0.5,0.9,2)*0.5, y=apply(d_pdist_j,2,max)[-1]*c(1,1,0.9,1,1), 
           label=c("G1","G2","G3","G4","G5"), colour=scales::hue_pal()(5), size=5, fontface="bold")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-35_pdist_j.emf", width = 7, height = 5)
g_pdist_j
dev.off()

msd_j <- matrix(c("theta","G1","G2","G3","G4","G5","mean",res_j$population_mean,"sd",res_j$population_sd),
                ncol = 6, byrow=T)

japa_cont <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/Vertical Scaling/japanese_content.csv", stringsAsFactors = F)



# math
res_m <- full_data_m %>% estip(fc=3, bg=3, ng=5, rm_list = c("e108", "e109"), min_a = 0)

res_m$para
res_m$SE

write.csv(res_m$para, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-12_itempara_m.csv", quote = F, row.names = F)
write.csv(res_m$SE, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-13_itempara_se_m.csv", quote = F, row.names = F)


eap_m <- full_data_m %>% estheta(param=res_m$para, mintheta = -4, maxtheta = 4, sigma = 3)

fit_m3 <- ifind(full_data_m, res_m$para, eap_m$res$EAP)
fit_m3$ggplot
fit_m3$X2
fit_m4 <- ifind3(full_data_m, res_m$para, eap_m$res$EAP)
write.csv(fit_m4, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-9-4_itemfit_m.csv", quote = F, row.names = F)


Q3_m2 <- Q3_stat(full_data_m, res_m$para, eap_m$res$EAP, abs_value = 0.3)
Q3_m2$df[Q3_m2$df$abs %>% isTRUE(),]

m <- res_m$para$a != 0 # key
res_m$para$level <- res_m$para$Item %>% str_extract("[a-z]{1,2}")
#ICC
d_icc_m <- apply(matrix(seq(-6,6,length.out = 101)), 1, ptheta, 
                 a=res_m$para$a[m], b=res_m$para$b[m], c=res_m$para$c[m], D=1.702) %>% t() %>% data.frame()
colnames(d_icc_m) <- res_m$para$Item[m]
d_icc_m$theta <- seq(-6,6,length.out = 101)
g_icc_m <- d_icc_m %>% tidyr::gather(key=Item, value=prob, -theta) %>% 
  ggplot(aes(x=theta, y=prob, group=Item, colour=Item))+
  geom_line()+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"), colour=TeX("Item"))
#IIF
d_iif_m <- apply(matrix(seq(-6,6,length.out = 101)), 1, iif, 
                 a=res_m$para$a[m], b=res_m$para$b[m], c=res_m$para$c[m], D=1.702) %>% t() %>% data.frame()
colnames(d_iif_m) <- res_m$para$Item[m]
d_iif_m$theta <- seq(-6,6,length.out = 101)
g_iif_m <- d_iif_m %>% tidyr::gather(key=Item, value=prob, -theta) %>% 
  ggplot(aes(x=theta, y=prob, group=Item, colour=Item))+
  geom_line()+#ylim(0,3.5)+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"), colour=TeX("Item"))

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-25_icc_,m.emf", width = 7, height = 5)
g_icc_m
dev.off()

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-26_iif_m.emf", width = 7, height = 5)
g_iif_m+ylim(0,3.8)
dev.off()

#scater
g_item_scat_m <- res_m$para[m,] %>% ggplot(aes(x=b, y=a, colour=Item))+
  geom_point()+xlim(-7,5)+ylim(0,2.5)+
  #scale_shape_manual(values=c(1,2,3,4,5,6,7))+
  labs(x="困難度", y="識別力", colour="Item")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-28_item_scatter_m.emf", width = 7, height = 5)
g_item_scat_m
dev.off()

#scatter 2
g_item_scat_m2 <- res_m$para[m,] %>% ggplot(aes(x=b, y=a, colour=level))+
  geom_point()+xlim(-7,5)+ylim(0,2.5)+
  geom_text_repel(aes(label=Item))+
  labs(x="困難度", y="識別力", colour="Item")+
  facet_wrap(~level, ncol = 1)+
  theme(legend.position = 'none')

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-33_item_scat_lv_m.emf", width = 7, height = 9)
g_item_scat_m2
dev.off()

# TIF
d_tif_m <- apply(matrix(seq(-6,6,length.out = 101)), 1, tif, 
                 a=res_m$para$a[m], b=res_m$para$b[m], c=res_m$para$c[m], D=1.702) %>% data.frame()
d_tif_m$theta <- seq(-6,6,length.out = 101)
g_tif_m <- d_tif_m %>%
  ggplot(aes(x=theta, y=.))+
  geom_line()+ylim(0,21)+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-30_tif_m.emf", width = 7, height = 5)
g_tif_m
dev.off()

# TIF of each levels
d_tif_lv_m <- data.frame(a=tif_auto(res_m$para[res_m$para$level=="a",], max = 10, min = -10),
                         b=tif_auto(res_m$para[res_m$para$level=="b",], max = 10, min = -10),
                         c=tif_auto(res_m$para[res_m$para$level=="c",], max = 10, min = -10),
                         d=tif_auto(res_m$para[res_m$para$level=="d",], max = 10, min = -10),
                         e=tif_auto(res_m$para[res_m$para$level=="e",], max = 10, min = -10),
                         f=tif_auto(res_m$para[res_m$para$level=="f",], max = 10, min = -10),
                         st=tif_auto(res_m$para[res_m$para$level=="st",], max = 10, min = -10),
                         theta=seq(-6,6,length.out = 301))
g_tif_lv_m <- d_tif_lv_m %>% tidyr::gather(key=levels, value=tif, -theta) %>% 
  ggplot(aes(x=theta, y=tif, colour=levels))+
  geom_line(aes(linetype=levels), size=1.5)+#ylim(0,6)+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))+
  annotate("text", x=c(-4,-2,2,3.7,-1,2.5,-1.2)*0.5, y=apply(d_tif_lv_m,2,max)[-8]*c(1,1,1,1,1,1,1), 
           label=c("a","b","c","d","e","f","st"), colour=scales::hue_pal()(7), size=5, fontface="bold")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-34_tif_lv_m.emf", width = 7, height = 5)
g_tif_lv_m
dev.off()

# population distribution
res_m$population_mean
d_pdist_m <- res_m$population_dist %>% data.frame() 
colnames(d_pdist_m) <- c("theta","G1","G2","G3","G4","G5")
g_pdist_m <- d_pdist_m %>% tidyr::gather(key=grade, value=Probability, -theta) %>% 
  ggplot(aes(x=theta, y=Probability, colour=grade))+
  geom_line()+
  annotate("text", x=c(-3,-2,-0.3,3.7,1)*0.5, y=apply(d_pdist_m,2,max)[-1]*c(1,1.1,1.05,0.7,1.1), 
           label=c("G1","G2","G3","G4","G5"), colour=scales::hue_pal()(5), size=5, fontface="bold")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-36_pdist_m.emf", width = 7, height = 5)
g_pdist_m
dev.off()

msd_m <- matrix(c("theta","G1","G2","G3","G4","G5","mean",res_m$population_mean,"sd",res_m$population_sd),
                ncol = 6, byrow=T)

math_cont <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/Vertical Scaling/math_content.csv", stringsAsFactors = F)


# output
# prior dist setting
# theta ~ Normal(0,1) # multigroup model
# b     ~ Normal(0,1.5^2)
# a     ~ logNormal(0, 0.25^2)

# Marginal Bayes
res_b_j <- full_data_j %>% estip(fc=3, bg=3, ng=5, min_a = 0, Bayes=1)
eap_b_j <- full_data_j %>% estheta(param = res_b_j$para)
res_b_m <- full_data_m %>% estip(fc=3, bg=3, ng=5, min_a = 0, Bayes=1, maxiter_em = 500) # 251回で収束
res_b_m <- full_data_m %>% estip(fc=3, bg=3, ng=5, min_a = 0, rm_list = c("a101"), Bayes=1, maxiter_em = 500)
eap_b_m <- full_data_m %>% estheta(param = res_b_m$para)

# parameter
res_b_j$para -> test
res_b_j$SE

res_b_m$para
res_b_m$SE

write.csv(res_b_j$para, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-17_itempara_b_j.csv", quote = F, row.names = F)
write.csv(res_b_j$SE, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-18_itempara_b_se_j.csv", quote = F, row.names = F)

write.csv(res_b_m$para, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-19_itempara_b_m.csv", quote = F, row.names = F)
write.csv(res_b_m$SE, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-20_itempara_b_se_m.csv", quote = F, row.names = F)


# fit
fit_b_j <- full_data_j%>% ifind(para=res_b_j$para, theta=eap_b_j$res$EAP, line_size = 0.5)
fit_b_j$ggplot
fit_b_j2 <- full_data_j%>% ifind3(para=res_b_j$para, theta=eap_b_j$res$EAP)

fit_b_m <- full_data_m%>% ifind(para=res_b_m$para, theta=eap_b_m$res$EAP, line_size = 0.5)
fit_b_m$ggplot
fit_b_m2 <- full_data_m%>% ifind3(para=res_b_m$para, theta=eap_b_m$res$EAP)

write.csv(fit_b_j2, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-20-1_itemfit_b_j.csv", quote = F, row.names = F)
write.csv(fit_b_m2, file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/table-for-M/4-20-2_itemfit_b_m.csv", quote = F, row.names = F)


devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-37_fit_graph_b_j.emf", width = 5, height = 9)
fit_b_j$ggplot+
  theme(axis.text = element_text(size=7),
        strip.text = element_text(size=7),
        strip.background = element_blank())
dev.off()

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-38_fit_graph_b_m.emf", width = 5, height = 9)
fit_b_m$ggplot+
  theme(axis.text = element_text(size=7),
        strip.text = element_text(size=7),
        strip.background = element_blank())
dev.off()

#japanese

#ICC
j <- res_b_j$para$a != 0 # key
res_b_j$para$level <- res_b_j$para$Item %>% str_extract("[a-z]{1,2}")

d_icc_b_j <- apply(matrix(seq(-6,6,length.out = 101)), 1, ptheta, 
                   a=res_b_j$para$a[j], b=res_b_j$para$b[j], c=res_b_j$para$c[j], D=1.702) %>% t() %>% data.frame()
colnames(d_icc_b_j) <- res_b_j$para$Item[j]
d_icc_b_j$theta <- seq(-6,6,length.out = 101)
g_icc_b_j <- d_icc_b_j %>% tidyr::gather(key=Item, value=prob, -theta) %>% 
  ggplot(aes(x=theta, y=prob, group=Item, colour=Item))+
  geom_line()+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"), colour=TeX("Item"))

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-39_icc_b_j.emf", width = 7, height = 5)
g_icc_b_j
dev.off()

#scatter 2
g_item_scat_b_j2 <- res_b_j$para[j,] %>% ggplot(aes(x=b, y=a, colour=level))+
  geom_point()+xlim(-5,5)+ylim(0,2)+
  geom_text_repel(aes(label=Item))+
  labs(x="困難度", y="識別力", colour="Item")+
  facet_wrap(~level, ncol = 1)+
  theme(legend.position = 'none')

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-40_item_scat_lv_b_j.emf", width = 7, height = 9)
g_item_scat_b_j2
dev.off()

# TIF of each levels
d_tif_lv_b_j <- data.frame(a=tif_auto(res_b_j$para[res_b_j$para$level=="a",]),
                           b=tif_auto(res_b_j$para[res_b_j$para$level=="b",]),
                           c=tif_auto(res_b_j$para[res_b_j$para$level=="c",]),
                           d=tif_auto(res_b_j$para[res_b_j$para$level=="d",]),
                           e=tif_auto(res_b_j$para[res_b_j$para$level=="e",]),
                           f=tif_auto(res_b_j$para[res_b_j$para$level=="f",]),
                           st=tif_auto(res_b_j$para[res_b_j$para$level=="st",]),
                           theta=seq(-6,6,length.out = 301))
g_tif_lv_b_j <- d_tif_lv_b_j %>% tidyr::gather(key=levels, value=tif, -theta) %>% 
  ggplot(aes(x=theta, y=tif, colour=levels))+
  geom_line(aes(linetype=levels), size=1.5)+ylim(0,6)+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))+
  annotate("text", x=c(-3.5,-2.5,-0.5,-1,-1.5,1,-1.5), y=apply(d_tif_lv_b_j,2,max)[-8]*c(1,1,1,0.9,1.1,1,0.8), 
           label=c("a","b","c","d","e","f","st"), colour=scales::hue_pal()(7), size=5, fontface="bold")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-41_tif_lv_j.emf", width = 7, height = 5)
g_tif_lv_b_j
dev.off()

# population distribution
res_b_j$population_mean
d_pdist_b_j <- res_b_j$population_dist %>% data.frame() 
colnames(d_pdist_b_j) <- c("theta","G1","G2","G3","G4","G5")
g_pdist_b_j <- d_pdist_b_j %>% tidyr::gather(key=grade, value=Probability, -theta) %>% 
  ggplot(aes(x=theta, y=Probability, colour=grade))+
  geom_line()+
  annotate("text", x=c(-1.5,0,0,2,1), y=apply(d_pdist_b_j,2,max)[-1]*c(1,1.05,0.95,0.4,1.1), 
           label=c("G1","G2","G3","G4","G5"), colour=scales::hue_pal()(5), size=5, fontface="bold")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-45_pdist_b_j.emf", width = 7, height = 5)
g_pdist_b_j
dev.off()

msd_b_j <- matrix(c("theta","G1","G2","G3","G4","G5","mean",res_b_j$population_mean,"sd",res_b_j$population_sd),
                ncol = 6, byrow=T)

#math

#ICC
m <- res_b_m$para$a != 0 # key
res_b_m$para$level <- res_b_m$para$Item %>% str_extract("[a-z]{1,2}")

d_icc_b_m <- apply(matrix(seq(-6,6,length.out = 101)), 1, ptheta, 
                 a=res_b_m$para$a[m], b=res_b_m$para$b[m], c=res_b_m$para$c[m], D=1.702) %>% t() %>% data.frame()
colnames(d_icc_b_m) <- res_b_m$para$Item[m]
d_icc_b_m$theta <- seq(-6,6,length.out = 101)
g_icc_b_m <- d_icc_b_m %>% tidyr::gather(key=Item, value=prob, -theta) %>% 
  ggplot(aes(x=theta, y=prob, group=Item, colour=Item))+
  geom_line()+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"), colour=TeX("Item"))

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-42_icc_b_m.emf", width = 7, height = 5)
g_icc_b_m
dev.off()

#scatter 2
g_item_scat_b_m2 <- res_b_m$para[m,] %>% ggplot(aes(x=b, y=a, colour=level))+
  geom_point()+xlim(-5,5)+ylim(0,2)+
  geom_text_repel(aes(label=Item))+
  labs(x="困難度", y="識別力", colour="Item")+
  facet_wrap(~level, ncol = 1)+
  theme(legend.position = 'none')

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-43_item_scat_lv_b_m.emf", width = 7, height = 9)
g_item_scat_b_m2
dev.off()

# TIF of each levels
d_tif_lv_b_m <- data.frame(a=tif_auto(res_b_m$para[res_b_m$para$level=="a",]),
                         b=tif_auto(res_b_m$para[res_b_m$para$level=="b",]),
                         c=tif_auto(res_b_m$para[res_b_m$para$level=="c",]),
                         d=tif_auto(res_b_m$para[res_b_m$para$level=="d",]),
                         e=tif_auto(res_b_m$para[res_b_m$para$level=="e",]),
                         f=tif_auto(res_b_m$para[res_b_m$para$level=="f",]),
                         st=tif_auto(res_b_m$para[res_b_m$para$level=="st",]),
                         theta=seq(-6,6,length.out = 301))
g_tif_lv_b_m <- d_tif_lv_b_m %>% tidyr::gather(key=levels, value=tif, -theta) %>% 
  ggplot(aes(x=theta, y=tif, colour=levels))+
  geom_line(aes(linetype=levels), size=1.5)+ylim(0,6)+
  labs(x=TeX("$\\theta$"), y=TeX("$I(\\theta)$"))+
  annotate("text", x=c(-2,-1.5,2,3,-0.8,5,-1), y=apply(d_tif_lv_b_m,2,max)[-8]*c(1,1,1,1,1,0.6,1), 
           label=c("a","b","c","d","e","f","st"), colour=scales::hue_pal()(7), size=5, fontface="bold")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-44_tif_lv_m.emf", width = 7, height = 5)
g_tif_lv_b_m
dev.off()

# population distribution
res_b_m$population_mean
d_pdist_b_m <- res_b_m$population_dist %>% data.frame() 
colnames(d_pdist_b_m) <- c("theta","G1","G2","G3","G4","G5")
g_pdist_b_m <- d_pdist_b_m %>% tidyr::gather(key=grade, value=Probability, -theta) %>% 
  ggplot(aes(x=theta, y=Probability, colour=grade))+
  geom_line()+
  annotate("text", x=c(-1.6,-1,0,2.3,1), y=apply(d_pdist_b_m,2,max)[-1]*c(1,1.05,1.1,0.4,1), 
           label=c("G1","G2","G3","G4","G5"), colour=scales::hue_pal()(5), size=5, fontface="bold")

devEMF::emf("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-46_pdist_b_m.emf", width = 7, height = 5)
g_pdist_b_m
dev.off()

msd_b_m <- matrix(c("theta","G1","G2","G3","G4","G5","mean",res_b_m$population_mean,"sd",res_b_m$population_sd),
                  ncol = 6, byrow=T)



# GIRT

res_GIRT <-  guIRT_fix(Uc=full_data_j, groupvar = "pop.num", idvar = "ID"
                       , ncat = rep(2,70), type = rep("B",70), 
          npointth=21, npointph=5, phmin=0, phmax=2, paramab=c(1,4), print=1 )


# oral examination

# icc
ggplot(data = data.frame(x=c(-4:4)), aes(x=x))+
  stat_function(fun=icc, args = list(a=0.5,b=0,D=1.702), aes(colour="a=0.5, b=0"), size=1.5)+
  stat_function(fun=icc, args = list(a=1.2,b=-2,D=1.702), aes(colour="a=1.2, b=-2"), size=1.5)+
  stat_function(fun=icc, args = list(a=2,b=1.5,D=1.702), aes(colour="a=2, b=1.5"), size=1.5)+
  labs(colour="discrimination\n & difficulty")+
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"))+ 
  theme_bw()

ggsave("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot_for_oral_exam/1-4_icc2pl.png", width = 7, height = 5)


d_pdist_b_m %>% tidyr::gather(key=grade, value=Probability, -theta) %>% 
  ggplot(aes(x=theta, y=Probability, colour=grade))+
  geom_line(size=c(rep(1.5, 31), rep(0.5, 31), rep(0.5, 31), rep(0.5, 31), rep(1.5, 31)))+
  annotate("text", x=c(-1.6,-1,0,2.3,1), y=apply(d_pdist_b_m,2,max)[-1]*c(1,1.05,1.1,0.4,1), 
           label=c("小4","小5","小6","中1","中2"), colour=scales::hue_pal()(5), size=5, fontface="bold")+
  theme_bw()+
  scale_color_hue(name="学年", labels=c("小4","小5","小6","中1","中2"))+
  labs(x=TeX("$\\theta$"))

ggsave("C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot_for_oral_exam/1-5_pdist.png", width = 8, height = 5)


