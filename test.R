one <- data.frame(matrix(c(1,2,3,4,5,6,7,8,9), 3, 3))

twe <- data.frame(matrix(c(9,8,7,6,5,4,3,2,1), 3, 3))

res <- one == twe

table(res)


# これなら，度数が0でも数え上げられる。
table(factor(res[,1], levels = c(TRUE, FALSE)))


# うまくいく関数がなさそうなので，自作関数を用意
Fun <- function(x){
  
  a <- sum(x == T)
  b <- sum(x == F)
  
  return(c(a,b))
}


apply(res, 2, Fun)



# 今度は，一致していない行数，列数を見つける。

grep(F,res[,1])

# これで良い。
apply(res, 2, grep, pattern = F)

