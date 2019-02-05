library(tidyverse)
library(irtfun2)


# gen simulation data set
th_1 <- rnorm(1000, 0, 1)
th_2 <- rnorm(1000, 0, 1)
th_3 <- rnorm(1000, 0, 1)

nofitems <- 70

true_a <- rlnorm(nofitems, 0, 0.4)
true_b <- rnorm(nofitems, 0, 1)

# LETTERS is built-in constatn, that has the 26 upper-case letters of the Roman alphabet
ID <- character(0)
for(t in 1:(nofitems/10)){
  LET <- LETTERS[t] %>% rep(10)
  NUM <- c(1:10) %>% formatC(width = 3, flag = 0) %>% as.matrix() 
  ID <- c(ID, paste0(LET, NUM))
}

t_para <- data.frame(Item = ID, a = true_a, b = true_b, c = rep(0, nofitems))

g1 <- ID[1:30]
g2 <- ID[21:50]
g3 <- ID[41:70]

sim_gen2 <- function(theta, para){
  a <- para$a
  b <- para$b
  c <- para$c
  ID <- para$Item
  dat <- sim_gen(theta = theta, a = a, b = b, c = c)
  #colnames(dat)[-1] <- ID
  return(dat)
}

dat1 <- t_para[ID %in% g1,] %>% sim_gen2(theta = th_1)
colnames(dat1)[-1] <- g1

dat2 <- t_para[ID %in% g2,] %>% sim_gen2(theta = th_2)
colnames(dat2)[-1] <- g2

dat3 <- t_para[ID %in% g3,] %>% sim_gen2(theta = th_3)
colnames(dat3)[-1] <- g3


res_1 <- dat1 %>% estip2(fc = 2)
res_2 <- dat2 %>% estip2(fc = 2)
res_3 <- dat3 %>% estip2(fc = 2)
res_1$para$Item <- res_1$para$Item %>% as.character()
res_2$para$Item <- res_2$para$Item %>% as.character()
res_3$para$Item <- res_3$para$Item %>% as.character()

#CEquating(T_para = res_1$para, F_para = res_2$para)

# 項目パラメタの平均と標準偏差を計算し，データフレームに納める
# 平均は項目困難度の平均を，標準偏差は外れ値を考慮して識別力幾何平均の逆数とする。

library(pracma)
pre_mean <- matrix(nrow = 3, ncol = 3)
pre_sd <- matrix(nrow = 3, ncol = 3)
# diag(pre_mean) <- c(res_1$para$b %>% mean(), res_2$para$b %>% mean(), res_2$para$b %>% mean())
# diag(pre_sd) <- c(res_1$para$a %>% geomean(), res_2$para$a %>% geomean(), res_2$para$a %>% geomean())

para_list <- list(res_1, res_2, res_3)
for(i in 1:3){ # of row
  row_a <- para_list[[i]]$para$a
  row_b <- para_list[[i]]$para$b
  row_item <- para_list[[i]]$para$Item
  for(j in 1:3){ # of column
    col_item <- para_list[[j]]$para$Item
    key <- row_item %in% col_item
    if(sum(key) == 0){
      pre_mean[i, j] <- NA
      pre_sd[i, j] <- NA
      next
    }
    target_a <- row_a[key]
    target_b <- row_b[key]
    pre_mean[i, j] <- target_b %>% mean()
    pre_sd[i, j] <- target_a %>% geomean()
    # pre_sd[i, j] <- target_b %>% sd()
  }
}

pre_mean <- matrix(c(68.4, 63.3, 0, 53.6, 60.2, 63.1, 55.3, 0, 0, 62.9, 57.1, 54.4, 76.5, 0, 58.9, 54.0), nrow=4)
pre_sd <- matrix(c(18.8, 3.4, 0, 5.4, 16.5, 3.5, 14.1, 0, 0, 3.6, 12.7, 6.8, 17.3, 0, 10.9, 6.1), nrow = 4)

# n of forms
# form <- c(1, 2, 3)
form <- c(1,2,3,4) # n of rater

# n of common subjects
n <- matrix(0, nrow = length(form), ncol = length(form))
for(i in form){
  row_item <- para_list[[i]]$para$Item
  for(j in form){
    col_item <- para_list[[j]]$para$Item
    n[i, j] <- row_item %in% col_item %>% sum()
  }
}

n <- matrix(c(20, 10, 0, 10, 10, 20, 10, 0, 0, 10, 20, 10, 10, 0, 10, 20), nrow = 4, ncol = 4)

# sigma = A coefficient
H1 <- matrix(0, nrow = length(form), ncol = length(form))

for(i in form){
  for(j in form){
    if(i == j){ # diag element
      for(k in form[!form == j]){
        if(is.na(pre_sd[i, k])) next
        H1[i, j] <- H1[i, j] + n[i, k] * pre_sd[i, k]^2
      }
    } else { # non diag
      H1[i, j] <- -n[i, j] * pre_sd[i, j] * pre_sd[j, i]
    }
    if(n[i, j] == 0){ # no common item between test two forms
      H1[i, j] <- 0
      next
    }
  }
}

inv_H1 <- eigen(H1, symmetric = F)
#inv_H1 <- pinv(H1)
A <- inv_H1$vectors[,length(form)] # objective vector
# A <- inv_H1[,length(form)]

# mean = K coefficient
H2 <- matrix(0, nrow = length(form), ncol = length(form))
for(i in form){
  for(j in form){
    if(n[i, j] == 0){ # no common item between test two forms
      H2[i, j] <- 0
      next
    }
    if(i == j){ # diag element
      for(k in form[!form == j]){
        if(is.na(pre_mean[i, k])) next
        H2[i, j] <- H2[i, j] + n[i, k] * pre_mean[i, k]
      }
    } else { # non diag
      H2[i, j] <- -n[i, j] * pre_mean[i, j] 
    }
  }
}


H3 <- matrix(0, nrow = length(form), ncol = length(form))
for(i in form){
  for(j in form){
    if(n[i, j] == 0){ # no common item between test two forms
      H3[i, j] <- 0
      next
    }
    if(i == j){ # diag element
      for(k in form[!form == j]){
        if(is.na(pre_mean[i, k])) next
        H3[i, j] <- H3[i, j] + n[i, k]
      }
    } else { # non diag
      H3[i, j] <- -n[i, j]
    }
  }
}

K <- -pinv(H3) %*% t(H2) %*% A
A
K
# res_1$para%>% dplyr::mutate(a = a/A[1]) %>% dplyr::mutate(b = A[1]*b + K[1])

# equate_manually <- function(para, A, K){
#   #para %>% dplyr::mutate(a = a/A) %>% dplyr::mutate(b = A*b + K)
#   para$a <- para$a/A
# }
