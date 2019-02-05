# Operation check of the program for IRT observed score

library(tidyverse
        )
# read ngaku data
ngaku16 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/02_文科委託研究/20170807 文科委託H29/2 板宮 得点分布の産出プログラム/ngakudata/ngaku16.csv", header=F)

# estimate item parameter
res_ngaku16 <- MML_EM_cpp(ngaku16,fc=2,ng=1,gc=0)

res_ngaku16$para

# estimate ability parameter
res_ngaku16$theta <- estheta(ngaku16,res_ngaku16$para[,c(2:4)],Gc=0,ITEMc=2, est="EAP")

# Check the estimater
#item parameter
res_ngaku16$para
# EAP
res_ngaku16$theta$res$EAP


# IRT observed score
obs_score <- score_dist(res_ngaku16$theta$res$EAP, res_ngaku16$para$a, res_ngaku16$para$b)

# rawscore
raw_score <- res_ngaku16$theta$res$SCORE

# true score
true_score <- tscore_dist(res_ngaku16$theta$res$EAP, res_ngaku16$para$a, res_ngaku16$para$b)

# freqency
obs_hist <- obs_score %>% table() %>% as.vector()
raw_hist <- raw_score %>% table() %>% as.vector()
true_hist <- true_score %>% table %>% as.vector()
true_hist <- c(0,true_hist,0,0) # 

# data.frame
hist_data <- data.frame(obs=obs_hist, raw=raw_hist, true=true_hist,score=0:25)

# tidy and ggplot
hist_data %>% tidyr::gather(key="type", value="freq", -score) %>% ggplot(aes(x=score,y=freq,colour=type))+
  geom_point()+
  geom_line()
