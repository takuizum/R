irf <- function(theta,a,b){
  1/(1+exp(-1.702*a*(theta-b)))
}

theta <- seq(-4,4,length.out = 301)

a <- c(1.0,1,1.5)
b <- c(-1.9,0,0)

df <- data.frame(theta=theta)
df$icc1 <- irf(theta,a[1],b[1])
df$icc2 <- irf(theta,a[2],b[2])
df$icc3 <- irf(theta,a[3],b[3])

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/2plm_example.emf",width=10,height = 5)
df %>% tidyr::gather(key=item, value=probability, -theta) %>% 
  ggplot(aes(x=theta, y=probability, group=item, linetype=item))+
  geom_line()
dev.off()


#######################
# Experiment 1
#######################

# ngaku data analysis
library(tidyverse)
library(irtfun2)

ngaku_all <- read_csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/12_R/R/ngakudata/ngaku-all.csv",
                      na=c("N"),col_names = F, col_types = str_flatten(rep("i",52)))

str(ngaku_all)
ngaku_all %>% purrr::map(unique)

# replace grade number to integer
ngaku_all$X2 %>% unique()

ngaku_all$X2[ngaku_all$X2==161] <- 1
ngaku_all$X2[ngaku_all$X2==182] <- 2
ngaku_all$X2[ngaku_all$X2==221] <- 3
ngaku_all$X2[ngaku_all$X2==222] <- 4
ngaku_all$X2[ngaku_all$X2==223] <- 5
ngaku_all$X2[ngaku_all$X2==224] <- 6

res_ngaku_all <- ngaku_all %>% estip(fc=3, eMLL = 1e-10)

i2pl <- function(theta,para,D){
  # 2pl only
  a <- para$a
  b <- para$b
  p <- 1/(1+exp(-D*a*(theta-b)))
  p
}
i2ti <- function(theta,para,D){
  a <- para$a
  b <- para$b
  p <- i2pl(theta,para,D)
  q <- 1-p
  D^2*a^2*p*q
}

#res_ngaku_all$para

# icc
theta_vec <- seq(-6,6,length.out = 301)
icc_df <- theta_vec %>% matrix() %>% apply(1,i2pl,para=res_ngaku_all$para,D=1.702) %>% 
  t() %>% as.data.frame()
icc_df$theta <- theta_vec

icc_g <- icc_df %>% tidyr::gather(key=item,value=probability,-theta) %>% 
  ggplot(aes(x=theta,y=probability,colour=item))+
  geom_line()+
  theme(legend.position = 'none')

# test info
tic_df <-theta_vec %>% matrix() %>% apply(1,i2ti,para=res_ngaku_all$para,D=1.702) %>% 
  t() %>% as.data.frame()
tic_df$theta <- theta_vec

tic_g <- tic_df %>% tidyr::gather(key=item,value=information,-theta) %>% 
  ggplot(aes(x=theta,y=information,colour=item))+
  geom_line()+
  theme(legend.position = 'none')

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/icc&iic.emf",width=10,height = 5)
gridExtra::grid.arrange(icc_g,tic_g,ncol=2)
dev.off()

# 並行テストに仕上がっているかどうかを，テスト情報量の観点から観察
key_1 <- ngaku_all[ngaku_all$X2==1,c(-1,-2)] %>% colSums() %>% is.na()
key_2 <- ngaku_all[ngaku_all$X2==2,c(-1,-2)] %>% colSums() %>% is.na()
key_5 <- ngaku_all[ngaku_all$X2==5,c(-1,-2)] %>% colSums() %>% is.na()
key_6 <- ngaku_all[ngaku_all$X2==6,c(-1,-2)] %>% colSums() %>% is.na()

tic_df_1 <-theta_vec %>% matrix() %>% apply(1,i2ti,para=res_ngaku_all$para[key_1,],D=1.702) %>% 
  t() %>% rowSums()
tic_df_2 <-theta_vec %>% matrix() %>% apply(1,i2ti,para=res_ngaku_all$para[key_2,],D=1.702) %>% 
  t() %>% rowSums()
tic_df_5 <-theta_vec %>% matrix() %>% apply(1,i2ti,para=res_ngaku_all$para[key_5,],D=1.702) %>% 
  t() %>% rowSums()
tic_df_6 <-theta_vec %>% matrix() %>% apply(1,i2ti,para=res_ngaku_all$para[key_6,],D=1.702) %>% 
  t() %>% rowSums()
tic_df <- data.frame(test1=tic_df_1,test2=tic_df_2,test3=tic_df_5,test4=tic_df_6,theta=theta_vec)

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/tic_4tests.emf",width=10,height = 5)
tic_df %>% tidyr::gather(key=test,value=information,-theta) %>% 
  ggplot(aes(x=theta,y=information,colour=test,linetype=test))+
  geom_line()
dev.off()
# IRT true score dist and observed score dist
# tibble format is not able to use in estheta finction.
# reload
ngaku_all <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/12_R/R/ngakudata/ngaku-all.csv", header = F,na.strings = "N")
res_theta <- ngaku_all %>% estheta(param=res_ngaku_all$para)
res_theta$res$GROUP %>% unique()
res_theta$res$GROUP %>% table()

# the comparison of distribution H16 and H18
# use all items
ngaku_ose_16 <- irtfun2::obscore_dist(theta=res_theta$res$EAP[res_theta$res$GROUP==161],
                                      a=res_ngaku_all$para$a, b=res_ngaku_all$para$b)
ngaku_ose_18 <- irtfun2::obscore_dist(theta=res_theta$res$EAP[res_theta$res$GROUP==182],
                                      a=res_ngaku_all$para$a, b=res_ngaku_all$para$b)

obs_df <- data.frame(score = 0:50, H16 = dist.f(ngaku_ose_16)$cum.pcnt, H18 = dist.f(ngaku_ose_18)$cum.pcnt)

obs_g <- obs_df %>% tidyr::gather(key=year, value=cumlative_percent, -score) %>% 
  ggplot(aes(x=score,y=cumlative_percent, colour=year))+
  geom_line()+
  geom_point()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_H16&H18.emf",width=10,height = 5)
obs_g
dev.off()

# the comparison obs VS true

ngaku_tse <- irtfun2::irt_ytx(res_ngaku_all$para[key_1,],res_ngaku_all$para[key_2,])

ngaku_ose_16 <- irtfun2::obscore_dist(theta=res_theta$res$EAP[res_theta$res$GROUP==161],
                                      a=res_ngaku_all$para$a[key_1], b=res_ngaku_all$para$b[key_1])
ngaku_ose_18 <- irtfun2::obscore_dist(theta=res_theta$res$EAP[res_theta$res$GROUP==182],
                                      a=res_ngaku_all$para$a[key_2], b=res_ngaku_all$para$b[key_2])

ngaku_ose <- irtfun2::epe(ngaku_ose_16,ngaku_ose_18)
#irtfun2::epe(ngaku_ose_18,ngaku_ose_16) # vice versa

se_df <- data.frame(X = ngaku_ose$X, irt_obs=ngaku_ose$eYx_U, irt_true=ngaku_tse$tau_Y)

se_g <- se_df %>% tidyr::gather(key=method,value=Y,-X) %>% 
  ggplot(aes(x=X,y=Y,colour=method))+
  geom_line()+
  geom_point()+
  labs(x="H16",y="H18")

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_obs&true.emf",width=7,height = 5)
se_g
dev.off()


#######################
# Experiment 2
#######################

ggplot(data.frame(a=0:4),aes(x=a))+
  stat_function(fun=dlnorm,args = list(mean=0,sd=0.5),colour="blue")+
  #stat_function(fun=dlnorm,args = list(mean=0.5,sd=1))+ # 平均を下げると，山がすこし平らになる。
  stat_function(fun=dlnorm,args = list(mean=-1,sd=0.5),colour="red")+ # 分散を小さくすると，山が正にゆがむ。
  stat_function(fun=dlnorm,args = list(mean=0.5,sd=0.2),colour="green")

# set.seed
set.seed(123)

theta_vec <- rnorm(5000)
a1 <- rlnorm(30,-1,0.5) # low
a2 <- rlnorm(30,0,0.5) # middle
a3 <- rlnorm(30,0.5,0.2) # high

b1 <- rnorm(30, mean=-1,sd=0.5) # low
b2 <- rnorm(30, mean=0, sd=1) # middle
b3 <- rnorm(30, mean=1, sd=0.5) # high
c <- rep(0,30) # asymptote

# condition 1
# stick slope middle
# low, middle, and high difficulty.

# set.seed
set.seed(123)
dat1 <- irtfun2::sim_gen(theta=theta_vec, a=a2,b=b1,c=c)
dat2 <- irtfun2::sim_gen(theta=theta_vec, a=a2,b=b2,c=c)
dat3 <- irtfun2::sim_gen(theta=theta_vec, a=a2,b=b3,c=c)

res_b1 <- dat1 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_b2 <- dat2 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_b3 <- dat3 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)

theta_b1 <- dat1 %>% estheta(param=res_b1$para,fc=2,gc=0)
theta_b2 <- dat2 %>% estheta(param=res_b2$para,fc=2,gc=0)
theta_b3 <- dat3 %>% estheta(param=res_b3$para,fc=2,gc=0)

b1_df <- data.frame(score=0:30,
                    raw=dist.f(dat1[,-1] %>% rowSums())$cddf,
                    obs_true = dist.f(obscore_dist(theta_vec,a2,b1))$cddf,
                    true_true = dist.f(tscore_dist(theta_vec,a2,b1))$cddf,
                    obs_est = dist.f(obscore_dist(theta_b1$res$EAP,res_b1$para$a,res_b1$para$b))$cddf,
                    true_est = dist.f(tscore_dist(theta_b1$res$EAP,res_b1$para$a,res_b1$para$b))$cddf %>% c(1)
                    )

b2_df <- data.frame(score=0:30,
                    raw=dist.f(dat2[,-1] %>% rowSums())$cddf,
                    obs_true = dist.f(obscore_dist(theta_vec,a2,b2))$cddf,
                    true_true = dist.f(tscore_dist(theta_vec,a2,b2))$cddf,
                    obs_est = dist.f(obscore_dist(theta_b2$res$EAP,res_b2$para$a,res_b2$para$b))$cddf,
                    true_est = dist.f(tscore_dist(theta_b2$res$EAP,res_b2$para$a,res_b2$para$b))$cddf %>% c(1,1)
)

b3_df <- data.frame(score=0:30,
                    raw=dist.f(dat3[,-1] %>% rowSums())$cddf,
                    obs_true = dist.f(obscore_dist(theta_vec,a2,b3))$cddf,
                    true_true = dist.f(tscore_dist(theta_vec,a2,b3))$cddf %>% c(1),
                    obs_est = dist.f(obscore_dist(theta_b3$res$EAP,res_b3$para$a,res_b3$para$b))$cddf,
                    true_est = dist.f(tscore_dist(theta_b3$res$EAP,res_b3$para$a,res_b3$para$b))$cddf %>% c(1,1)
)

b1_g <- b1_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-M, location-L" )+
  theme(plot.title=element_text(size=13,hjust=0.5),
        legend.position = 'none')

b2_g <- b2_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-M, location-M" )+
  theme(plot.title=element_text(size=13,hjust=0.5),
        legend.position = 'none')

b3_g <- b3_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-M, location-H" )+
  theme(plot.title=element_text(size=13,hjust=0.5),
        legend.position = 'none')



# condition 2
# stick slope low
# low, middle, and high difficulty

# set.seed
set.seed(123)
dat4 <- irtfun2::sim_gen(theta=theta_vec, a=a1,b=b1,c=c)
dat5 <- irtfun2::sim_gen(theta=theta_vec, a=a1,b=b2,c=c)
dat6 <- irtfun2::sim_gen(theta=theta_vec, a=a1,b=b3,c=c)

res_a1 <- dat4 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_a2 <- dat5 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_a3 <- dat6 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)

theta_a1 <- dat4 %>% estheta(param=res_a1$para,fc=2,gc=0)
theta_a2 <- dat5 %>% estheta(param=res_a2$para,fc=2,gc=0)
theta_a3 <- dat6 %>% estheta(param=res_a3$para,fc=2,gc=0)

a1_df <- data.frame(score=0:30,
                    raw=dist.f(dat4[,-1] %>% rowSums())$cddf,
                    obs_true = dist.f(obscore_dist(theta_vec,a1,b1))$cddf,
                    true_true = dist.f(tscore_dist(theta_vec,a1,b1))$cddf %>% c(1,1,1),
                    obs_est = dist.f(obscore_dist(theta_a1$res$EAP,res_a1$para$a,res_a1$para$b))$cddf,
                    true_est = dist.f(tscore_dist(theta_a1$res$EAP,res_a1$para$a,res_a1$para$b))$cddf %>% c(rep(1,4))
)

a2_df <- data.frame(score=0:30,
                    raw=dist.f(dat5[,-1] %>% rowSums())$cddf %>% c(1),
                    obs_true = dist.f(obscore_dist(theta_vec,a1,b2))$cddf %>% c(1),
                    true_true = dist.f(tscore_dist(theta_vec,a1,b2))$cddf %>% c(rep(1,4)),
                    obs_est = dist.f(obscore_dist(theta_a2$res$EAP,res_a2$para$a,res_a2$para$b))$cddf %>% c(1),
                    true_est = dist.f(tscore_dist(theta_a2$res$EAP,res_a2$para$a,res_a2$para$b))$cddf %>% c(rep(1,6))
)

a3_df <- data.frame(score=0:30,
                    raw=dist.f(dat6[,-1] %>% rowSums())$cddf %>% c(1,1),
                    obs_true = dist.f(obscore_dist(theta_vec,a1,b3))$cddf %>% c(rep(1,3)),
                    true_true = dist.f(tscore_dist(theta_vec,a1,b3))$cddf %>% c(rep(1,6)),
                    obs_est = dist.f(obscore_dist(theta_a3$res$EAP,res_a3$para$a,res_a3$para$b))$cddf %>% c(rep(1,4)),
                    true_est = dist.f(tscore_dist(theta_a3$res$EAP,res_a3$para$a,res_a3$para$b))$cddf %>% c(rep(1,7))
)

a1_g <- a1_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-L, location-L" )+
  theme(plot.title=element_text(size=13,hjust=0.5),
        legend.position = 'none')

a2_g <- a2_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-L, location-M" )+
  theme(plot.title=element_text(size=13,hjust=0.5),
        legend.position = 'none')

a3_g <- a3_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-L, location-H" )+
  theme(plot.title=element_text(size=13,hjust=0.5),
        legend.position = 'none')



# condition 3
# stick slope high
# low, middle, and high difficulty

# set.seed
set.seed(123)
dat7 <- irtfun2::sim_gen(theta=theta_vec, a=a3,b=b1,c=c)
dat8 <- irtfun2::sim_gen(theta=theta_vec, a=a3,b=b2,c=c)
dat9 <- irtfun2::sim_gen(theta=theta_vec, a=a3,b=b3,c=c)

res_c1 <- dat7 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_c2 <- dat8 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_c3 <- dat9 %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)

theta_c1 <- dat7 %>% estheta(param=res_c1$para,fc=2,gc=0)
theta_c2 <- dat8 %>% estheta(param=res_c2$para,fc=2,gc=0)
theta_c3 <- dat9 %>% estheta(param=res_c3$para,fc=2,gc=0)

c1_df <- data.frame(score=0:30,
                    raw=dist.f(dat7[,-1] %>% rowSums())$cddf,
                    obs_true = dist.f(obscore_dist(theta_vec,a3,b1))$cddf,
                    true_true = dist.f(tscore_dist(theta_vec,a3,b1))$cddf,
                    obs_est = dist.f(obscore_dist(theta_c1$res$EAP,res_c1$para$a,res_c1$para$b))$cddf,
                    true_est = dist.f(tscore_dist(theta_c1$res$EAP,res_c1$para$a,res_c1$para$b))$cddf
)

c2_df <- data.frame(score=0:30,
                    raw=dist.f(dat8[,-1] %>% rowSums())$cddf,
                    obs_true = dist.f(obscore_dist(theta_vec,a3,b2))$cddf,
                    true_true = dist.f(tscore_dist(theta_vec,a3,b2))$cddf,
                    obs_est = dist.f(obscore_dist(theta_c2$res$EAP,res_c2$para$a,res_c2$para$b))$cddf,
                    true_est = dist.f(tscore_dist(theta_c2$res$EAP,res_c2$para$a,res_c2$para$b))$cddf
)

c3_df <- data.frame(score=0:30,
                    raw=dist.f(dat9[,-1] %>% rowSums())$cddf,
                    obs_true = dist.f(obscore_dist(theta_vec,a3,b3))$cddf,
                    true_true = dist.f(tscore_dist(theta_vec,a3,b3))$cddf,
                    obs_est = dist.f(obscore_dist(theta_c3$res$EAP,res_c3$para$a,res_c3$para$b))$cddf,
                    true_est = dist.f(tscore_dist(theta_c3$res$EAP,res_c3$para$a,res_c3$para$b))$cddf %>% c(1)
)

c1_g <- c1_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-H, location-L" )+
  theme(plot.title=element_text(size=13,hjust=0.5),
        legend.position = 'none')

c2_g <- c2_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-H, location-M" )+
  theme(plot.title=element_text(size=13,hjust=0.5),
        legend.position = 'none')

c3_g <- c3_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,shape=method,linetype=method))+
  geom_line()+
  geom_point()+
  labs(title="slope-H, location-H" )+
  theme(plot.title=element_text(size=13,hjust=0.5))


emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison-true&obs-recover.emf",width=7,height = 5)
gridExtra::grid.arrange(a1_g,a2_g,a3_g,
                        b1_g,b2_g,b3_g,
                        c1_g,c2_g,c3_g,ncol=3)
dev.off()

# グラフが潰れてしまって分かりにくいので，個別に表示したグラフも保存
emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aLbL.emf",width=7,height = 5)
a1_g + theme(legend.position = "right")
dev.off()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aLbM.emf",width=7,height = 5)
a2_g + theme(legend.position = "right")
dev.off()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aLbH.emf",width=7,height = 5)
a3_g + theme(legend.position = "right")
dev.off()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aMbL.emf",width=7,height = 5)
b1_g + theme(legend.position = "right")
dev.off()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aMbM.emf",width=7,height = 5)
b2_g + theme(legend.position = "right")
dev.off()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aMbH.emf",width=7,height = 5)
b3_g + theme(legend.position = "right")
dev.off()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aHbL.emf",width=7,height = 5)
c1_g + theme(legend.position = "right")
dev.off()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aHbM.emf",width=7,height = 5)
c2_g + theme(legend.position = "right")
dev.off()

emf(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_aHbH.emf",width=7,height = 5)
c3_g + theme(legend.position = "right")
dev.off()



# facet data frame
comp_all_df <- data.frame(score=tidyr::gather(data=a1_df,key=method,value="a-L_b-L",-score)$score)
comp_all_df$method = tidyr::gather(data=a1_df,key=method,value="a-L_b-L",-score)$method
comp_all_df$"a-L_b-L" = tidyr::gather(data=a1_df,key=method,value="a-L_b-L",-score)$"a-L_b-L"
comp_all_df$"a-L_b-M" = tidyr::gather(data=a2_df,key=method,value="a-L_b-M",-score)$"a-L_b-M"
comp_all_df$"a-L_b-H" = tidyr::gather(data=a3_df,key=method,value="a-L_b-H",-score)$"a-L_b-H"
comp_all_df$"a-M_b-L" = tidyr::gather(data=b1_df,key=method,value="a-M_b-L",-score)$"a-M_b-L"
comp_all_df$"a-M_b-M" = tidyr::gather(data=b2_df,key=method,value="a-M_b-M",-score)$"a-M_b-M"
comp_all_df$"a-M_b-H" = tidyr::gather(data=b3_df,key=method,value="a-M_b-H",-score)$"a-M_b-H"
comp_all_df$"a-H_b-L" = tidyr::gather(data=c1_df,key=method,value="a-H_b-L",-score)$"a-H_b-L"
comp_all_df$"a-H_b-M" = tidyr::gather(data=c2_df,key=method,value="a-H_b-M",-score)$"a-H_b-M"
comp_all_df$"a-H_b-H" = tidyr::gather(data=c3_df,key=method,value="a-H_b-H",-score)$"a-H_b-H"

comp_all_df <- comp_all_df %>% tidyr::gather(key=cond, value=cddf, -score, -method)

comp_all_df$cond_f <- factor(comp_all_df$cond, levels = c("a-L_b-L","a-L_b-M","a-L_b-H","a-M_b-L","a-M_b-M","a-M_b-H",
                                                          "a-H_b-L","a-H_b-M","a-H_b-H"))
comp_all_df$method <- factor(comp_all_df$method, levels = c("raw", "true_true", "obs_true", "true_est", "obs_est"))


comp_all_g <- comp_all_df%>% 
  ggplot(aes(x=score,y=cddf,colour=method,linetype=method))+
  geom_line()+
  facet_wrap(~cond_f, nrow=3, labeller = label_parsed)

emf(file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_all.emf",width=7,height = 5)
comp_all_g
dev.off()

  
#######################
# Experiment 3
#######################

# Experiment2 で生成したデータを一部流用。
# 平均的な学力分布，困難度，識別力で，受検者の条件を変更させた場合
# set.seed
set.seed(123)
theta_vecA <- rnorm(100)
theta_vecB <- rnorm(500)
theta_vecC <- rnorm(1000)
theta_vecD <- theta_vec
theta_vecE <- rnorm(10000)

# 使用する項目パラメタはngakuデータとも比較できるように，middle slope & low location

#dat1 <- irtfun2::sim_gen(theta=theta_vec, a=a2,b=b1,c=c)

dat_A <- irtfun2::sim_gen(theta=theta_vecA,a=a2,b=b1,c=c)
dat_B <- irtfun2::sim_gen(theta=theta_vecB,a=a2,b=b1,c=c)
dat_C <- irtfun2::sim_gen(theta=theta_vecC,a=a2,b=b1,c=c)
dat_D <- irtfun2::sim_gen(theta=theta_vecD,a=a2,b=b1,c=c)
dat_E <- irtfun2::sim_gen(theta=theta_vecE,a=a2,b=b1,c=c)

res_A <- dat_A %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_B <- dat_B %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_C <- dat_C %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_D <- dat_D %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)
res_E <- dat_E %>% estip(fc=2,gc=0,e_ell = 0,min_a = 0.01, EM_dist = 0)

# 推定に失敗した項目の識別力を0と代入し，推定しないようにする。
res_A$para$a[10] <- 0

theta_A <- dat_A %>% estheta(param=res_A$para,fc=2,gc=0)
theta_B <- dat_B %>% estheta(param=res_B$para,fc=2,gc=0)
theta_C <- dat_C %>% estheta(param=res_C$para,fc=2,gc=0)
theta_D <- dat_D %>% estheta(param=res_D$para,fc=2,gc=0)
theta_E <- dat_E %>% estheta(param=res_E$para,fc=2,gc=0)

A_df <- data.frame(score=0:30,
                    raw=dist.f(dat_A[,-1] %>% rowSums())$cddf,
                    obs_true = dist.f(obscore_dist(theta_vecA,a2,b1))$cddf,
                    true_true = dist.f(tscore_dist(theta_vecA,a2,b1))$cddf %>% c(1),
                    obs_est = dist.f(obscore_dist(theta_A$res$EAP,res_A$para$a,res_A$para$b))$cddf,
                    true_est = dist.f(tscore_dist(theta_A$res$EAP,res_A$para$a,res_A$para$b))$cddf %>% c(1,1)
)

B_df <- data.frame(score=0:30,
                   raw=dist.f(dat_B[,-1] %>% rowSums())$cddf,
                   obs_true = dist.f(obscore_dist(theta_vecB,a2,b1))$cddf,
                   true_true = dist.f(tscore_dist(theta_vecB,a2,b1))$cddf,
                   obs_est = dist.f(obscore_dist(theta_B$res$EAP,res_B$para$a,res_B$para$b))$cddf,
                   true_est = dist.f(tscore_dist(theta_B$res$EAP,res_B$para$a,res_B$para$b))$cddf %>% c(1)
)

C_df <- data.frame(score=0:30,
                   raw=dist.f(dat_C[,-1] %>% rowSums())$cddf,
                   obs_true = dist.f(obscore_dist(theta_vecC,a2,b1))$cddf,
                   true_true = dist.f(tscore_dist(theta_vecC,a2,b1))$cddf,
                   obs_est = dist.f(obscore_dist(theta_C$res$EAP,res_C$para$a,res_C$para$b))$cddf,
                   true_est = dist.f(tscore_dist(theta_C$res$EAP,res_C$para$a,res_C$para$b))$cddf %>% c(1)
)

D_df <- data.frame(score=0:30,
                   raw=dist.f(dat_D[,-1] %>% rowSums())$cddf,
                   obs_true = dist.f(obscore_dist(theta_vecD,a2,b1))$cddf,
                   true_true = dist.f(tscore_dist(theta_vecD,a2,b1))$cddf,
                   obs_est = dist.f(obscore_dist(theta_D$res$EAP,res_D$para$a,res_D$para$b))$cddf,
                   true_est = dist.f(tscore_dist(theta_D$res$EAP,res_D$para$a,res_D$para$b))$cddf %>% c(1)
)

E_df <- data.frame(score=0:30,
                   raw=dist.f(dat_E[,-1] %>% rowSums())$cddf,
                   obs_true = dist.f(obscore_dist(theta_vecE,a2,b1))$cddf,
                   true_true = dist.f(tscore_dist(theta_vecE,a2,b1))$cddf,
                   obs_est = dist.f(obscore_dist(theta_E$res$EAP,res_E$para$a,res_E$para$b))$cddf,
                   true_est = dist.f(tscore_dist(theta_E$res$EAP,res_E$para$a,res_E$para$b))$cddf %>% c(1)
)

A_g <- A_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,linetype=method))+
  geom_line()+
  labs(title="N-100" )+
  theme(plot.title=element_text(size=13,hjust=0.5))

B_g <- B_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,linetype=method))+
  geom_line()+
  labs(title="N-500" )+
  theme(plot.title=element_text(size=13,hjust=0.5))

C_g <- C_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,linetype=method))+
  geom_line()+
  labs(title="N-1000" )+
  theme(plot.title=element_text(size=13,hjust=0.5))

D_g <- D_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,linetype=method))+
  geom_line()+
  labs(title="N-5000" )+
  theme(plot.title=element_text(size=13,hjust=0.5))

E_g <- E_df %>% tidyr::gather(key=method,value=cddf,-score) %>% 
  ggplot(aes(x=score,y=cddf,colour=method,linetype=method))+
  geom_line()+
  labs(title="N-10000" )+
  theme(plot.title=element_text(size=13,hjust=0.5))


comp_sub_df <- A_df %>% tidyr::gather(key=method,value="N-100",-score)
comp_sub_df$"N-500" <- tidyr::gather(data=B_df,key=method,value="N-500",-score)$"N-500"
comp_sub_df$"N-1000" <- tidyr::gather(data=C_df,key=method,value="N-1000",-score)$"N-1000"
comp_sub_df$"N-5000" <- tidyr::gather(data=D_df,key=method,value="N-5000",-score)$"N-5000"
comp_sub_df$"N-10000" <- tidyr::gather(data=E_df,key=method,value="N-10000",-score)$"N-10000"

comp_sub_df <- comp_sub_df %>% tidyr::gather(key=cond,value=cddf,-score,-method)

comp_sub_df$cond <- factor(comp_sub_df$cond,levels = c("N-100","N-500","N-1000","N-5000","N-10000"))

comp_sub_df$method <- factor(comp_sub_df$method, levels = c("raw", "true_true", "obs_true", "true_est", "obs_est"))

comp_sub_g <- comp_sub_df %>% ggplot(aes(x=score,y=cddf,colour=method,linetype=method))+
  geom_line()+
  facet_wrap(~cond)

emf(file="C:/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/20181011_plot-for-report/comparison_sub.emf",width=7,height = 5)
comp_sub_g
dev.off()
