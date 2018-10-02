# MCMC


# MCMC1

N <- 1e4 # サンプリング回数

MCMC1_R <- function(N){
  
  a <- numeric(N)　#　代入用ベクトル
  b <- rep(NA,N)
  x <- 0
  
  model <- function(x){
    1/(1 + x^2)
  }
  
  p <- model(x)
  
  accept <- 0
  
  for(i in 1:N){
    
    # prior distribution
    # this distribution must be symmetric.
    y <- runif(1,x-5, x+5)
    q <- model(y)
    
    if(runif(1) < q/p){
      
      x <- y
      p <- q
      accept <- accept + 1
      
    } else {
      
      b[i] <- y
      
    }
    
    a[i] <- x
    #cat("iteration",i,"\r")
    
  }
  
  res <- data.frame("accept" = a, "reject" = b)
  
  return(res)
}


accept/N

hist(a, breaks = 100)

system.time(
  MCMC1_R(1e6)
)


# Rcppで記述したファイルを読み込み，こちらでの実行速度を確認

sourceCpp(file = "C:/Users/GenericSSC/OneDrive/Documents/12_R/Rcpp/20180707_model_cpp.cpp")

sourceCpp(file = "C:/Users/sep10_z1vk0al/OneDrive/Documents/12_R/Rcpp/20180707_MCMC1_cpp.cpp")


system.time(
  MCMC1_cpp(1e5)
)

res <- MCMC1_cpp(1e6)

hist(res$accept, breaks = 10000, xlim = c(-10,10))


# cppの方が5倍は早い。
# MCMC2



