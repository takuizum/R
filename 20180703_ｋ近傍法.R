head(iris)

# check data mean
summary(iris)


# k-最近傍法
dat_learn <- sample_frac(tbl = iris, size = 0.6)



dat <-split(iris,list(c("learn","predict")))


model <- kknn(Species ~ ., dat$learn, dat$predict)

summary(model)
fitted(model)


# k-近傍法
library(class)
res <- knn(train = dat$learn[,-5], test = dat$predict[,-5], cl = dat$learn[,5] ,k = 3)

table(res)

sum(res == dat$predict[,5])/length(res)



pairs(iris[,1:4], pch=as.character(iris[,5]), col = c(2,3,4))
