library(extraDistr)

plot(x = seq(0.0001, 4, length.out = 100), y = dbetapr(seq(0.0001, 4, length.out = 100), 1, 3))

smpl <- rbetapr(1000000, 1, 3)
M <- mean(smpl)
V <- var(smpl)
beta <- M*(M+1)/V + 2
alpha <- M*(beta - 1)
rm(list = ls())

bppar <- function(prob, node){
  if(round(sum(prob), digits = 5) != 1) prob <- prob/sum(prob)
  # cat("Theoritical expectation and variance is", , )
  M <- node %*% prob ; print(as.vector(M))
  V <- (node - c(M))^2 %*% prob ; print(as.vector(V))
  beta <- M*(M+1)/V + 2
  alpha <- M*(beta - 1)
  c("shape1" = alpha, "shape2" = beta)
}


bppar(node = seq(0.0001, 4, length.out = 30), prob = dbetapr(seq(0.0001, 4, length.out = 30), 1, 3))
E <- 1/(3 - 1)
V <- 1*(3 + 1 -1)/((3 - 2)*(3 - 1)^2)
