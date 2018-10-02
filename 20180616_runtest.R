sample(c(0,1,2,NA), 10, replace = T, prob = c(0.15, 0.6, 0.15, 0.1))

num <- as.data.frame(rep(50,100))

# apply(num, 1, sample, x = c(0,1,2,NA), replace = T, prob = c(0.15, 0.6, 0.15, 0.1))


# create sample data set

dat <- t(
  apply(num, 1, sample, x = c(0,1,2,NA), replace = T, prob = c(0.15, 0.6, 0.15, 0.1))
)

str(dat)

dat1 <- dat

dat2 <- dat
# matrix
system.time(
  dat[dat == 2] <- 0
)



# data.table

dat <- as.data.table(dat)

system.time(
  dat[dat == 2] <- 0
)


# dplyr
dat <- tbl_df(dat)

system.time(
  dat[dat == 2] <- 0
)



# 行列で実行するのが，一番早い。



# replacement test

replace(dat[1,], which(dat[1, ] %in% 2), 0)


# これでも置き換えはできるがキャラクター型になるのが，厄介。
gsub(2, 0, dat[1,])


# これがめちゃくちゃ簡単だった。
dat[dat == 2] <- 0
