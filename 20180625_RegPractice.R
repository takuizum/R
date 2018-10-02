setwd("C:/Users/sep10_z1vk0al/OneDrive/Documents/TEMP")


dat1 <- read.csv("JyugyouHyouka.csv", header = T)

dat1[,-1]

apply(dat1[,-1], 2, mean)
apply(dat1[,-1], 2, sd)


mr <- lm(y~x1+x2+x3+x4, data = dat1)
summary(mr)


library(car)

vif(mr)


mrZ <- data.frame(scale(dat1[,-1]))
mrZ <- lm(y~x1+x2+x3+x4, data = mrZ)
summary(mrZ, digits=3)
vif(mrZ)


library(MASS)


# 変数増加法による変数選択
mrAIC <- stepAIC(lm(y~1, data = dat1), direction = "forward", scope = list(upper = ~x1+x2+x3+x4))

# 変数減少法による変数選択
mrAIC <- stepAIC(lm(y~x1+x2+x3+x4, data = dat1), direction = "backward", scope = list(upper = ~1))

# Step wise
mrAIC <- stepAIC(lm(y~1, data = dat1), direction = "both", scope = list(upper = ~x1+x2+x3+x4, lowwer = ~1))


# 変数選択した内容を反映させて，再度重回帰分析


dat2 <- read.csv("20180618testdata.csv", na.strings = "N")
dat2 <- na.omit(dat2)
rownames(dat2[,1])

因子番号 <- c("insi",rep(1,5), rep(4,8), 3, 3, rep(2,7), rep(3,4), rep(2,5), rep(3,5),2,rep(1,7), 2)

dat2 <- rbind(dat2,因子番号)

id <- dat2[,1]
res <- numeric(0)

for(i in 1:4){
  
  key <- 因子番号 == i
  
  sub <- dat2[,key]
  
  res1 <- rowSums(sub)
  
  res <- cbind(res, res1)

}

 colnames(res) <- c("Enjo", "Gakusyu", "Syousatsu", "Rakkan")

 apply(res, 2, mean)
 apply(res, 2, sd) 

 res <- as.data.frame(scale(res))
 
 model <- Rakkan ~ Enjo + Gakusyu + Syousatsu
 
 Reg <- lm(model, data = res)
 
 # Step wise
 mrAIC <- stepAIC(Reg, direction = "both", scope = list(upper = ~Enjo+Gakusyu+Syousatsu, lowwer = ~1))
 
 
 
 # final
 
 Reg <- lm(Rakkan~ Enjo+Gakusyu, data = res)
 
