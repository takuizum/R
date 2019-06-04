# batch treatment

library(mirt)
library(tidyverse)
library(EstCRM)

category <- c(2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 30, 40, 50)
subject <- c(1000, 5000, 10000) 
mctime <- 50

# evaluation methods
# 1. paremeter RMSE(theta, discrimination, mean of difficulty, logistic transformation of observed score)
# 2. AIC
# 

# using functions
# item sampling
grm_para_gen <- function(category_vector, b_dist = "norm", a_dist = "lnorm", 
                         args_b_dist = list(0, 1), args_a_dist = list(0, 1), 
                         item_name = NULL, index = "Item"){
  assign("rand_dist1", get(paste0("r", b_dist), mode = "function"))
  assign("rand_dist2", get(paste0("r", a_dist), mode = "function"))
  if(a_dist == "beta"){
    res <- purrr::map(as.list(category_vector - 1), 
                      function(x){c(rand_dist2(1, args_a_dist[[1]], args_a_dist[[2]])*1.5+0.4, 
                                    sort(rand_dist1(n = x, args_b_dist[[1]], args_b_dist[[2]])))})
  } else {
    res <- purrr::map(as.list(category_vector - 1), 
                      function(x){c(rand_dist2(1, args_a_dist[[1]], args_a_dist[[2]]), 
                                    sort(rand_dist1(n = x, args_b_dist[[1]], args_b_dist[[2]])))})
  }
  if(is.null(item_name)){
    item_name <- formatC(c(1:length(category_vector)), width = 3, flag = 0)
    item_name <- as.vector(apply(matrix(index, ncol = 1), 1, paste0, item_name))
  }
  names(res) <- item_name
  return(res)
}

# Probability density function
pcrm <- function(theta, y, xi, a, b, c){ # c = 1/alpha
  z <- log(y/(xi-y))
  bz <- b + c * z
  Const <- xi/(1/c*(y*xi-y^2))
  dnorm(theta, mean = bz, sd = 1/a) * Const
}
lpcrm <- function(theta, y, xi, a, b, c){ # c = 1/alpha
  z <- log(y/(xi-y))
  bz <- b + c * z
  Const1 <- a*xi*c
  Const2 <- sqrt(2*pi)*(y*xi-y^2)
  log(Const1/Const2)-a^2/2*(theta - bz)^2
}
Lik_crm <- function(theta, y, xi, lambda, width){
  a <- lambda[,1]
  b <- lambda[,2]
  c <- lambda[,3]
  prod(pcrm(theta, y, xi, a, b, c)*width, na.rm = T)
}
Lik_crm2 <- function(theta, y, xi, lambda, const){ # all_pattern must be matrix with 1 column.
  a <- lambda[,1]
  b <- lambda[,2]
  c <- lambda[,3]
  prod(pcrm(theta, y, xi, a, b, c) / const, na.rm = T)
}
mle <- function(z, para){
  a <- para[,1]
  b <- para[,2]
  c <- para[,3]
  bz <- b + c * z
  sum(a^2 * bz, na.rm = T) / sum(a^2, na.rm = T)
}
logit_adj_aic <- function(x, max_cat , min_cat, delta, transformation = "GRM"){ # logit変換による補正項だけを求める
  # x <- data_long
  if(transformation == "GRM"){
    X <- x$resp
    X <- x$resp * (max_cat[x$item] - min_cat[x$item] + delta + delta)
    X <- X + min_cat[x$item] - delta
    XI <- max_cat[x$item] + delta
  } else if(transformation == "W&Z"){
    X <- x$resp * (max_cat[x$item] - min_cat[x$item])
    if (length(which(X == 0.01)) != 0) {
      X[which(X == 0.01)] <-  0
    }
    if (length(which(X == (max_cat[x$item] - min_cat[x$item] - 0.01))) != 0) {
      X[which(X == (max_cat[x$item] - min_cat[x$item] - 0.01))] <- 
        (max_cat[x$item] - min_cat[x$item])[which(X == (max_cat[x$item] - min_cat[x$item] - 0.01))]
    }
    X <- X + 1
    XI <- max_cat[x$item] + delta
  } else {
    NULL
  }
  -2*sum(log( 1/(XI - X) )) # adjustment term
}

estCRip2 <- function(x, fc = 3, IDc = 0, Gc = 0, Ntheta = 31, delta = 1, eEM = 0.001, eM = 0.001, eDIST = 1e-3, 
                     theta_range = c(-4, 4), transformation = "GRM", print = 1, maxiter_em = 100, 
                     est_dist = FALSE){
  # use data
  data <- as_tibble(x[,fc:ncol(x)])
  max_cat <- unlist(map(data, max))
  min_cat <- unlist(map(data, min))
  zdata <- data
  if(transformation == "GRM"){
    # suitable for GRM type transformation
    delta <- 1
    for(j in 1:ncol(data)){
      data[, j] <- data[, j] - min_cat[j] + delta # 最小のカテゴリは0（未観測）
      data[, j] <- data[, j] / (max_cat[j] - min_cat[j] + delta + delta) # 最大のカテゴリは1（未観測）
      zdata[, j] <- log(data[, j] / (1 - data[, j])) # logit transformation
    }
  } else if(transformation == "W&Z"){
    # # EstCRM type transformation
    for (i in 1:ncol(data)) {
      data[, i] <- data[, i] - min_cat[i]
      if (length(which(data[, i] == (max_cat[i] - min_cat[i]))) != 0) {
        data[which(data[, i] == (max_cat[i] - min_cat[i])), i] = (max_cat[i] - min_cat[i]) - 0.01
      }
      if (length(which(data[, i] == 0)) != 0) {
        data[which(data[, i] == 0), i] = 0.01
      }
      data[, i] <- data[, i] / (max_cat[i] - min_cat[i])
      zdata[, i] <- log(data[, i] / (1 - data[, i]))
    }
  }
  
  J <- ncol(data) # of items
  N <- nrow(data) # of subjects
  # prior weight
  M <- Ntheta
  Xq <- seq(theta_range[1], theta_range[2], length.out = Ntheta)
  Aq <- dnorm(Xq)/sum(dnorm(Xq))
  # initial value
  t0 <- cbind(1/sqrt(abs(unlist(map_df(as_tibble(zdata), var))-1)),
              -unlist(map_df(as_tibble(zdata), mean)),
              rep(1, J))
  
  #
  data <- as.matrix(data)
  zdata <- as.matrix(zdata)
  Lim <- matrix(0, ncol = M, nrow = nrow(data))
  colnames(data) <- c(1:J)
  data_long <- data %>% as_tibble %>% gather(key = item, value = resp) %>% map_df(as.numeric)
  
  tEM <- 0
  convEM <- TRUE
  while(convEM){
    tEM <- tEM + 1
    if(print > 1) cat(tEM, "time iteration.\n")
    #//////////////////////
    # E-step
    #//////////////////////
    mu1 <- numeric(nrow(data))
    sigma1 <-sqrt( 1/(sum(t0[,1]^2, na.rm = T) + 1^(-2)))
    for(i in 1:nrow(data)){
      # In future, this code will be rewrited using apply function
      mu1[i] <- sigma1^2 * sum(t0[,1]^2 * (t0[,2] + t0[,3] * zdata[i,]), na.rm = T)
    }
    term2 <- 0
    for(i in 1:N){
      term2[i] <- sum(t0[,1]^2 * (t0[,2] + t0[,3] * zdata[i,] - mu1[i])^2 + sigma1^2, na.rm = T)
    }
    ell <- N * sum(log(t0[,1]) + log(t0[,3])) - 0.5 * sum(term2) - ((N * J/2) * log(2 * pi)) # except constant term
    if(print <= 1) cat(sprintf("Expected Log Likelihood %.5f \b",ell), "\r")
    if(print > 1) cat(sprintf("Expected Log Likelihood %.5f \b",ell), "\n")    #//////////////////////
    # M-step 
    #//////////////////////
    t1 <- t0
    mm <- mean(mu1, na.rm = T)
    vm <- sum(mu1^2)/N - mm^2
    for(j in 1:J){
      z <- zdata[,j]
      mz <- mean(z, na.rm = T)
      vz <- sum(z^2, na.rm = T) / N - mz^2 
      czm <- sum(z * mu1, na.rm = T) / N - mz*mm
      t1[j,3] <- (vm + sigma1^2) / czm # c
      t1[j,2] <- mm - t1[j,3] * mz # b
      t1[j,1] <- (t1[j,3]^2*vz - vm - sigma1^2)^(-0.5)
    } # end of M step
    
    # scale transformation
    if(print > 1) cat("Empirical population distribution: mu =", mm, "sigma =", sqrt(vm), "\n")
    # equating coefficient
    if(tEM > 10000){ # Because empirical distribution estimated by initial value is very unstable, the caribration was not done
      A <- 1 / sqrt(vm)
      K <- 0 - A*mm
      # a
      t1[1,] <- t1[1,]/A
      t1[2,] <- t1[2,]*A + K
      t1[3,] <- t1[3,]*A # c is inverse of alpha
    }
    # EM cycle convergence
    if( max(abs(t1-t0)) < eEM && tEM != 1) convEM <- FALSE
    
    t0 <- t1
  }
  cat("\nConverged!\n")
  # AIC calculate
  # normalization constant
  all_pattern <- data %>% as_tibble %>% map(unique) %>% map(sort)
  const <- matrix(0, nrow = J, ncol = M)
  for(j in 1:J){
    xj <- matrix(all_pattern[[j]])
    for(m in 1:M){
      const[j,m] <- sum(apply(xj, 1, pcrm, theta = Xq[m], xi = 1, a = t0[j,1], b = t0[j,2], c = t0[j,3]))
    }
  }
  # Normalized Likelihood
  for(i in 1:nrow(data)){
    for(q in 1:M){
      Lim[i,q] <- Lik_crm2(Xq[q], data[i,], xi = 1, t0, const[,m])
    }
  }
  # loglikelihood
  lnL <- sum(log(Lim %*% Aq))
  # return list object
  colnames(t0) <- c("a", "b", "gamma")
  adj_term <- logit_adj_aic(data_long, max_cat = max_cat, min_cat = min_cat, delta = delta, transformation = transformation)
  # result
  list(para = as_tibble(t0),AIC = -2 * ell + 2*J*3 ,AIC_adj = -2 * ell + 2*J*3 + adj_term,
       AIC_norm = -2*lnL + 2*J*3, AIC_norm_adj = -2*lnL + adj_term + 2*J*3,lnL = ell, lnL_norm = lnL, adj_term = adj_term)
}
esthetaCR <- function(x, para, fc, IDc = NULL, transformation = "GRM", Ntheta = 31, theta_range = c(-4, 4)){
  # use data
  data <- as_tibble(x[,fc:ncol(x)])
  if(is.null(IDc)){
    ID <- 1:nrow(x)
  } else {
    ID <- x[,IDc]
  }
  # data(SelfEff)
  # data <- SelfEff
  max_cat <- unlist(map(data, max))
  min_cat <- unlist(map(data, min))
  zdata <- data
  
  if(transformation == "GRM"){
    # suitable for GRM type transformation
    delta <- 1
    for(j in 1:ncol(data)){
      data[, j] <- data[, j] - min_cat[j] + delta # 最小のカテゴリは0（未観測）
      data[, j] <- data[, j] / (max_cat[j] + delta) # 最大のカテゴリは1（未観測）
      zdata[, j] <- log(data[, j] / (1 - data[, j])) # logit transformation
    }
  } else if(transformation == "W&Z"){
    # # EstCRM type transformation
    for (i in 1:ncol(data)) {
      data[, i] <- data[, i] - min_cat[i]
      if (length(which(data[, i] == (max_cat[i] - min_cat[i]))) != 0) {
        data[which(data[, i] == (max_cat[i] - min_cat[i])), i] = (max_cat[i] - min_cat[i]) - 0.01
      }
      if (length(which(data[, i] == 0)) != 0) {
        data[which(data[, i] == 0), i] = 0.01
      }
      data[, i] <- data[, i] / (max_cat[i] - min_cat[i])
      zdata[, i] <- log(data[, i] / (1 - data[, i]))
    }
  }
  
  J <- ncol(data)
  # prior weight
  M <- Ntheta
  Xq <- seq(theta_range[1], theta_range[2], length.out = Ntheta)
  Aq <- dnorm(Xq)/sum(dnorm(Xq))
  width <- (theta_range[2]-theta_range[1])/(M-1)
  # initial value
  t0 <- as.matrix(para)
  data <- as.matrix(data)
  colnames(data) <- c(1:J)
  Lim <- matrix(0, ncol = M, nrow = nrow(data))
  Gim <- matrix(0, ncol = M, nrow = nrow(data))
  colnames(Gim) <- c(Xq)
  for(i in 1:nrow(data)){
    for(q in 1:M){
      Lim[i,q] <- Lik_crm(Xq[q], data[i,], xi = 1, t0, width)
    }
  }
  for(i in 1:nrow(data)){
    Gim[i, ] <- Lim[i,]*Aq / sum(Lim[i,]*Aq)
  }
  
  # EAP estimate
  cat("Calculating EAP...\n")
  eap <- Gim %*% Xq
  dev <- matrix(Xq, ncol = M, nrow = nrow(data), byrow = T) - matrix(rep(eap, M), ncol = M)
  se_eap <- sqrt(rowSums(dev^2 * Gim))
  
  # MLE estimate
  cat("Calculating MLE...\n")
  mle <- apply(zdata, 1, mle, para = para)
  se_mle <- rep(sqrt(1/(sum(para[,1]^2))), nrow(data))
  tibble(ID = ID, MLE = mle, SE_MLE = se_mle, EAP = as.vector(eap), SE_EAP = se_eap)
}



# mirtmodel
mod <- mirt.model("F1 = 1-10")

# nested list
subject_condition <- list_along(subject)
names(subject_condition) <- paste("subject", subject, sep = "_")

sink("/Users/takuizum/OneDrive/Documents/console.txt")
for(n in subject){
  category_condition <- list_along(category)
  names(category_condition) <- paste("category", category, sep = "_")
  for(k in category){
    cat("Estimating",k,"category data...\n\n\n")
    mctime_trial <- list_along(1:mctime)
    names(mctime_trial) <- paste("mctime", 1:mctime, sep = "_")
    mct <- 0
    while(mct != 50){ # stop when the counter reachs 50 
      # data generation routine----------------------------
      mct <- mct + 1
      cat(mct, "time simulation...\n")
      flag <- TRUE
      t <- 0
      while(flag){
        t <- t + 1
        seed <- round(runif(1)*1000000000)
        set.seed(seed)
        theta_sample <- rnorm(n, 0, 1) # scale metric
        item_cate <- rep(k, 10)
        item_sample <- grm_para_gen(item_cate, a_dist = "beta", b_dist = "norm",
                                    args_b_dist = list("mean" = 0, "sd" = 2), 
                                    args_a_dist = list("shape1" = 2,"shape2" = 2))
        a <- map_df(item_sample, function(x)x[1]) %>% t
        d <- map_df(item_sample, function(x) -x[-1]*x[1]) %>% t
        if(k == 2)
          dat <- mirt::simdata(a =  a, d = d, Theta = matrix(theta_sample), itemtype = "2PL", mins = 1)
        else 
          dat <- mirt::simdata(a =  a, d = d, Theta = matrix(theta_sample), itemtype = "graded", mins = 1)
        
        min_n <- dat %>% as_tibble %>% map(table) %>% map(min) %>% unlist %>% min
        cat_n <- dat %>% as_tibble %>% map(table) %>% map(length) %>% unlist %>% min
        
        cat(t, " times random sampling. The minimum N is", min_n,"and the number of categories is", cat_n,".\r")
        if(min_n > 1 && min_n != 0 && cat_n == k){
          cat("\nGood sample data has been found!\n")
          flag <- FALSE
        } 
      }
      # parameter estimation step
      if(k == 2)
        fit <- mirt(dat, mod, itemtype = "2PL") 
      else
        fit <- mirt(dat, mod, itemtype = "graded")
      grmpar <- coef(fit, IRTpars = T)
      grmtheta_mle <- fscores(fit, method = "ML", ) %>% as.vector()
      grmtheta_eap <- fscores(fit, method = "EAP") %>% as.vector()
      # CRM estimate
      res <- estCRip2(dat, fc = 1, transformation = "GRM")
      crmpar <- res$para
      res_theta <- esthetaCR(dat, crmpar, fc = 1)
      crmtheta_mle <- res_theta$MLE
      crmtheta_eap <- res_theta$EAP
      # EstCRM
      cat("Runnning EstCRM...")
      estcrm <- try(suppressWarnings(EstCRMitem(dat %>% as.data.frame, max.item = rep(k, 10), min.item = rep(1, 10))), silent = T)
      if(class(estcrm) == "try-error"){
        next
        mct <- mct - 1
      }
        # estcrm$LL
      crmtheta_mle2 <- try(EstCRMperson(dat %>% as.data.frame, ipar = estcrm$param, max.item = rep(k, 10), min.item = rep(1, 10)))
      # AIC
      aic_crm <- res$AIC
      aic_crm_adj <- res$AIC_adj
      aic_grm1 <- fit@Fit$AIC
      
      # par(mfrow=c(1,2)) 
      # plot(grmtheta_eap, crmtheta_eap, xlim = c(-4, 4), ylim = c(-4, 4))
      # plot(grmtheta_mle, crmtheta_mle, xlim = c(-4, 4), ylim = c(-4, 4))
      
      cat("\nAIC:  
     CRM    |  CRM(adj)  |  EstCRM   | GRM(mirt)   |
  ",aic_crm,"|  ",aic_crm_adj, "|", -2*estcrm$LL+2*3*15," |", aic_grm1,"   |\n\n\n")
      #----------------------------------------------
      tmp <- try(list(crmpar = crmpar, 
                  grmpar = grmpar, 
                  crmtheta = tibble(EAP = crmtheta_eap, MLE = crmtheta_mle, MLE2 = crmtheta_mle2$thetas[,2]), 
                  grmtheta = tibble(EAP = grmtheta_eap, MLE = grmtheta_mle), 
                  aic = c(crm = aic_crm, crm_adj = aic_crm_adj, EstCRM = -2*estcrm$LL+2*3*15, grm1 = aic_grm1),
                  seed = seed,
                  truetheta = theta_sample,
                  truepar = item_sample
                  ))
      
      mctime_trial[[paste("mctime", mct, sep = "_")]] <- tmp
    }
    category_condition[[paste("category", k, sep = "_")]] <- mctime_trial
  }
  subject_condition[[paste("subject", n, sep = "_")]] <- category_condition
}
sink()
# save
# save(list = c("subject_condition", "category", "subject", "mctime"), file = "grm_sim4.Rdata")
save(list = c("subject_condition", "category", "subject", "mctime"), file = "/Users/takuizum/OneDrive/Documents/grm_sim4.Rdata")
