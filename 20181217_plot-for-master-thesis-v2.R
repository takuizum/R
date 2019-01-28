# 修論シミュレーションプロット用
library(tidyverse)  # ggplot & tidyr data
library(RColorBrewer)  # usefull colour chart
library(devEMF)  # for output emf

cols <- brewer.pal(6, "Paired")

sim_list <- list("sim_all6-3")
#sep10_z1vk0al

for(sim in sim_list){
  setwd(paste0("/Users/sep10_z1vk0al/OneDrive/Documents/03_JSPS/Vertical Scaling/",sim))
  # 項目パラメタ，母集団平均，標準偏差ファイルの読み込み
  cat("reading data\r")
  cc <- 0
  calr_dicc_pdist <- numeric(9)
  sep_dicc_pdist <- numeric(9)
  con_dicc_pdist <- numeric(9)
  calr_dicc <- numeric(9)
  sep_dicc <- numeric(9)
  con_dicc <- numeric(9)
  itemsize <- numeric(9)
  subjectsize <- numeric(9)
  condition <- c("A1","A2","A3","B1","B2","B3","C1","C2","C3")
  for(sampleN in c(400,1000,10000)){
    for(itemN in c(15,30,60)){
      cc <- cc + 1
      # item parameter error
      cone_para <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_ParaError_con.csv"))
      sepe_para <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_ParaError_sep.csv"))
      calre_para <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_ParaError_calr.csv"))
      # pdist
      con_mean <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_pdist_mean_con.csv"))
      con_sd <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_pdist_sd_con.csv"))
      sep_mean <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_pdist_mean_sep.csv"))
      sep_sd <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_pdist_sd_sep.csv"))
      calr_mean <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_pdist_mean_calr.csv"))
      calr_sd <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_pdist_sd_calr.csv"))
      # DiCC(pdist)
      con_dicc_temp <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_DICC_con.csv"))
      sep_dicc_temp <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_DICC_sep.csv"))
      calr_dicc_temp <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_DICC_calr.csv"))
      # assign
      assign(sprintf("con_para_%d_%d",sampleN,itemN), cone_para)
      assign(sprintf("sep_para_%d_%d",sampleN,itemN), sepe_para)
      assign(sprintf("calr_para_%d_%d",sampleN,itemN), calre_para)
      
      assign(sprintf("con_mean_%d_%d",sampleN,itemN), con_mean)
      assign(sprintf("sep_mean_%d_%d",sampleN,itemN), sep_mean)
      assign(sprintf("calr_mean_%d_%d",sampleN,itemN), calr_mean)
      
      assign(sprintf("con_sd_%d_%d",sampleN,itemN), con_sd)
      assign(sprintf("sep_sd_%d_%d",sampleN,itemN), sep_sd)
      assign(sprintf("calr_sd_%d_%d",sampleN,itemN), calr_sd)
      
      calr_dicc_pdist[cc] <- mean(calr_dicc_temp$V1)
      sep_dicc_pdist[cc]  <- mean(sep_dicc_temp$V1)
      con_dicc_pdist[cc]  <- mean(con_dicc_temp$V1)
      
      # DICC
      calr_t <- sep_t <- con_t <- 0 # initialization
      for(t in 1:100){
        true <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_simN",t,"para_true.csv"), header = T)
        calr <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_simN",t,"para_calr.csv"), header = T)
        sep  <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_simN",t,"para_sep.csv"), header = T)
        con  <- read.csv(paste0("sampleN",sampleN,"_itemN",itemN,"_simN",t,"para_con.csv"), header = T)
        # theta interval and point
        N <- 31
        theta <- seq(-4,4, length.out = N) # Arai & Mayekawa (2011) では-3~3の区間だった。
        # 全項目，全分点の平均確率を計算。
        calr_t <- calr_t + mean(mapply(FUN = function(a1,a2,b1,b2,theta) mean(abs(1/(1+exp(-1.702*a2*(theta-b2))) - 1/(1+exp(-1.702*a1*(theta-b1))))), 
                                       a1=true$V3, a2=calr$p1, b1=true$V4, b2=calr$p2, MoreArgs = list(theta=theta), SIMPLIFY = T))
        
        sep_t  <- sep_t + mean(mapply(FUN = function(a1,a2,b1,b2,theta) mean(abs(1/(1+exp(-1.702*a2*(theta-b2))) - 1/(1+exp(-1.702*a1*(theta-b1))))), 
                                      a1=true$V3, a2=sep$V3, b1=true$V4, b2=sep$V4, MoreArgs = list(theta=theta), SIMPLIFY = T))
        
        con_t  <- con_t + mean(mapply(FUN = function(a1,a2,b1,b2,theta) mean(abs(1/(1+exp(-1.702*a2*(theta-b2))) - 1/(1+exp(-1.702*a1*(theta-b1))))), 
                                      a1=true$V3, a2=con$a, b1=true$V4, b2=con$b, MoreArgs = list(theta=theta), SIMPLIFY = T))
      } # end of t
      calr_dicc[cc] <- calr_t/100
      sep_dicc[cc]  <- sep_t/100
      con_dicc[cc]  <- con_t/100
      itemsize[cc]  <- itemN
      subjectsize[cc] <- sampleN
    }
  }
  
  rm(true,calr,sep,con,calr_t,sep_t,con_t,con_dicc_temp,sep_dicc_temp,calr_dicc_temp)
  rm(cone_para,sepe_para,calre_para,con_mean,sep_mean,calr_mean,con_sd,sep_sd,calr_sd)
  
  cat("create data.frame\r")
  # para
  # a
  RMSEa_sep <- data.frame(i400j15 = sep_para_400_15$RMSE.a, i400s30 = sep_para_400_30$RMSE.a, i400s60 = sep_para_400_60$RMSE.a, 
                          i1000j15 = sep_para_1000_15$RMSE.a, i1000j30 = sep_para_1000_30$RMSE.a, i1000j60 = sep_para_1000_60$RMSE.a,
                          i10000j15 = sep_para_10000_15$RMSE.a, i10000j30 = sep_para_10000_30$RMSE.a, i10000j60 = sep_para_10000_60$RMSE.a)
  
  RMSEa_calr <- data.frame(i400j15 = calr_para_400_15$RMSE.a, i400s30 = calr_para_400_30$RMSE.a, i400s60 = calr_para_400_60$RMSE.a, 
                           i1000j15 = calr_para_1000_15$RMSE.a, i1000j30 = calr_para_1000_30$RMSE.a, i1000j60 = calr_para_1000_60$RMSE.a,
                           i10000j15 = calr_para_10000_15$RMSE.a, i10000j30 = calr_para_10000_30$RMSE.a, i10000j60 = calr_para_10000_60$RMSE.a)
  
  RMSEa_con <- data.frame(i400j15 = con_para_400_15$RMSE.a, i400s30 = con_para_400_30$RMSE.a, i400s60 = con_para_400_60$RMSE.a, 
                          i1000j15 = con_para_1000_15$RMSE.a, i1000j30 = con_para_1000_30$RMSE.a, i1000j60 = con_para_1000_60$RMSE.a,
                          i10000j15 = con_para_10000_15$RMSE.a, i10000j30 = con_para_10000_30$RMSE.a, i10000j60 = con_para_10000_60$RMSE.a)
  
  #b
  RMSEb_sep <- data.frame(i400j15 = sep_para_400_15$RMSE.b, i400s30 = sep_para_400_30$RMSE.b, i400s60 = sep_para_400_60$RMSE.b, 
                          i1000j15 = sep_para_1000_15$RMSE.b, i1000j30 = sep_para_1000_30$RMSE.b, i1000j60 = sep_para_1000_60$RMSE.b,
                          i10000j15 = sep_para_10000_15$RMSE.b, i10000j30 = sep_para_10000_30$RMSE.b, i10000j60 = sep_para_10000_60$RMSE.b)
  
  RMSEb_calr <- data.frame(i400j15 = calr_para_400_15$RMSE.b, i400s30 = calr_para_400_30$RMSE.b, i400s60 = calr_para_400_60$RMSE.b, 
                           i1000j15 = calr_para_1000_15$RMSE.b, i1000j30 = calr_para_1000_30$RMSE.b, i1000j60 = calr_para_1000_60$RMSE.b,
                           i10000j15 = calr_para_10000_15$RMSE.b, i10000j30 = calr_para_10000_30$RMSE.b, i10000j60 = calr_para_10000_60$RMSE.b)
  
  RMSEb_con <- data.frame(i400j15 = con_para_400_15$RMSE.b, i400s30 = con_para_400_30$RMSE.b, i400s60 = con_para_400_60$RMSE.b, 
                          i1000j15 = con_para_1000_15$RMSE.b, i1000j30 = con_para_1000_30$RMSE.b, i1000j60 = con_para_1000_60$RMSE.b,
                          i10000j15 = con_para_10000_15$RMSE.b, i10000j30 = con_para_10000_30$RMSE.b, i10000j60 = con_para_10000_60$RMSE.b)
  # mean
  #sep
  sep_mean_g1 <- data.frame(i400j15=sep_mean_400_15$V1, i400s30 = sep_mean_400_30$V1, i400s60 = sep_mean_400_60$V1, 
                            i1000j15 = sep_mean_1000_15$V1, i1000j30 = sep_mean_1000_30$V1, i1000j60 = sep_mean_1000_60$V1,
                            i10000j15 = sep_mean_10000_15$V1, i10000j30 = sep_mean_10000_30$V1, i10000j60 = sep_mean_10000_60$V1)
  sep_mean_g2 <- data.frame(i400j15=sep_mean_400_15$V2, i400s30 = sep_mean_400_30$V2, i400s60 = sep_mean_400_60$V2, 
                            i1000j15 = sep_mean_1000_15$V2, i1000j30 = sep_mean_1000_30$V2, i1000j60 = sep_mean_1000_60$V2,
                            i10000j15 = sep_mean_10000_15$V2, i10000j30 = sep_mean_10000_30$V2, i10000j60 = sep_mean_10000_60$V2)
  sep_mean_g3 <- data.frame(i400j15=sep_mean_400_15$V3, i400s30 = sep_mean_400_30$V3, i400s60 = sep_mean_400_60$V3, 
                            i1000j15 = sep_mean_1000_15$V3, i1000j30 = sep_mean_1000_30$V3, i1000j60 = sep_mean_1000_60$V3,
                            i10000j15 = sep_mean_10000_15$V3, i10000j30 = sep_mean_10000_30$V3, i10000j60 = sep_mean_10000_60$V3)
  sep_mean_g4 <- data.frame(i400j15=sep_mean_400_15$V4, i400s30 = sep_mean_400_30$V4, i400s60 = sep_mean_400_60$V4, 
                            i1000j15 = sep_mean_1000_15$V4, i1000j30 = sep_mean_1000_30$V4, i1000j60 = sep_mean_1000_60$V4,
                            i10000j15 = sep_mean_10000_15$V4, i10000j30 = sep_mean_10000_30$V4, i10000j60 = sep_mean_10000_60$V4)
  sep_mean_g5 <- data.frame(i400j15=sep_mean_400_15$V5, i400s30 = sep_mean_400_30$V5, i400s60 = sep_mean_400_60$V5, 
                            i1000j15 = sep_mean_1000_15$V5, i1000j30 = sep_mean_1000_30$V5, i1000j60 = sep_mean_1000_60$V5,
                            i10000j15 = sep_mean_10000_15$V5, i10000j30 = sep_mean_10000_30$V5, i10000j60 = sep_mean_10000_60$V5)
  sep_mean_all <- dplyr::bind_rows(sep_mean_g1,sep_mean_g2,sep_mean_g3,sep_mean_g4,sep_mean_g5)
  #calr
  calr_mean_g1 <- data.frame(i400j15=calr_mean_400_15$V1, i400s30 = calr_mean_400_30$V1, i400s60 = calr_mean_400_60$V1, 
                             i1000j15 = calr_mean_1000_15$V1, i1000j30 = calr_mean_1000_30$V1, i1000j60 = calr_mean_1000_60$V1,
                             i10000j15 = calr_mean_10000_15$V1, i10000j30 = calr_mean_10000_30$V1, i10000j60 = calr_mean_10000_60$V1)
  calr_mean_g2 <- data.frame(i400j15=calr_mean_400_15$V2, i400s30 = calr_mean_400_30$V2, i400s60 = calr_mean_400_60$V2, 
                             i1000j15 = calr_mean_1000_15$V2, i1000j30 = calr_mean_1000_30$V2, i1000j60 = calr_mean_1000_60$V2,
                             i10000j15 = calr_mean_10000_15$V2, i10000j30 = calr_mean_10000_30$V2, i10000j60 = calr_mean_10000_60$V2)
  calr_mean_g3 <- data.frame(i400j15=calr_mean_400_15$V3, i400s30 = calr_mean_400_30$V3, i400s60 = calr_mean_400_60$V3, 
                             i1000j15 = calr_mean_1000_15$V3, i1000j30 = calr_mean_1000_30$V3, i1000j60 = calr_mean_1000_60$V3,
                             i10000j15 = calr_mean_10000_15$V3, i10000j30 = calr_mean_10000_30$V3, i10000j60 = calr_mean_10000_60$V3)
  calr_mean_g4 <- data.frame(i400j15=calr_mean_400_15$V4, i400s30 = calr_mean_400_30$V4, i400s60 = calr_mean_400_60$V4, 
                             i1000j15 = calr_mean_1000_15$V4, i1000j30 = calr_mean_1000_30$V4, i1000j60 = calr_mean_1000_60$V4,
                             i10000j15 = calr_mean_10000_15$V4, i10000j30 = calr_mean_10000_30$V4, i10000j60 = calr_mean_10000_60$V4)
  calr_mean_g5 <- data.frame(i400j15=calr_mean_400_15$V5, i400s30 = calr_mean_400_30$V5, i400s60 = calr_mean_400_60$V5, 
                             i1000j15 = calr_mean_1000_15$V5, i1000j30 = calr_mean_1000_30$V5, i1000j60 = calr_mean_1000_60$V5,
                             i10000j15 = calr_mean_10000_15$V5, i10000j30 = calr_mean_10000_30$V5, i10000j60 = calr_mean_10000_60$V5)
  calr_mean_all <- dplyr::bind_rows(calr_mean_g1,calr_mean_g2,calr_mean_g3,calr_mean_g4,calr_mean_g5)
  #con
  con_mean_g1 <- data.frame(i400j15=con_mean_400_15$V1, i400s30 = con_mean_400_30$V1, i400s60 = con_mean_400_60$V1, 
                            i1000j15 = con_mean_1000_15$V1, i1000j30 = con_mean_1000_30$V1, i1000j60 = con_mean_1000_60$V1,
                            i10000j15 = con_mean_10000_15$V1, i10000j30 = con_mean_10000_30$V1, i10000j60 = con_mean_10000_60$V1)
  con_mean_g2 <- data.frame(i400j15=con_mean_400_15$V2, i400s30 = con_mean_400_30$V2, i400s60 = con_mean_400_60$V2, 
                            i1000j15 = con_mean_1000_15$V2, i1000j30 = con_mean_1000_30$V2, i1000j60 = con_mean_1000_60$V2,
                            i10000j15 = con_mean_10000_15$V2, i10000j30 = con_mean_10000_30$V2, i10000j60 = con_mean_10000_60$V2)
  con_mean_g3 <- data.frame(i400j15=con_mean_400_15$V3, i400s30 = con_mean_400_30$V3, i400s60 = con_mean_400_60$V3, 
                            i1000j15 = con_mean_1000_15$V3, i1000j30 = con_mean_1000_30$V3, i1000j60 = con_mean_1000_60$V3,
                            i10000j15 = con_mean_10000_15$V3, i10000j30 = con_mean_10000_30$V3, i10000j60 = con_mean_10000_60$V3)
  con_mean_g4 <- data.frame(i400j15=con_mean_400_15$V4, i400s30 = con_mean_400_30$V4, i400s60 = con_mean_400_60$V4, 
                            i1000j15 = con_mean_1000_15$V4, i1000j30 = con_mean_1000_30$V4, i1000j60 = con_mean_1000_60$V4,
                            i10000j15 = con_mean_10000_15$V4, i10000j30 = con_mean_10000_30$V4, i10000j60 = con_mean_10000_60$V4)
  con_mean_g5 <- data.frame(i400j15=con_mean_400_15$V5, i400s30 = con_mean_400_30$V5, i400s60 = con_mean_400_60$V5, 
                            i1000j15 = con_mean_1000_15$V5, i1000j30 = con_mean_1000_30$V5, i1000j60 = con_mean_1000_60$V5,
                            i10000j15 = con_mean_10000_15$V5, i10000j30 = con_mean_10000_30$V5, i10000j60 = con_mean_10000_60$V5)
  con_mean_all <- dplyr::bind_rows(con_mean_g1,con_mean_g2,con_mean_g3,con_mean_g4,con_mean_g5)
  
  # sd
  #sep
  sep_sd_g1 <- data.frame(i400j15=sep_sd_400_15$V1, i400s30 = sep_sd_400_30$V1, i400s60 = sep_sd_400_60$V1, 
                          i1000j15 = sep_sd_1000_15$V1, i1000j30 = sep_sd_1000_30$V1, i1000j60 = sep_sd_1000_60$V1,
                          i10000j15 = sep_sd_10000_15$V1, i10000j30 = sep_sd_10000_30$V1, i10000j60 = sep_sd_10000_60$V1)
  sep_sd_g2 <- data.frame(i400j15=sep_sd_400_15$V2, i400s30 = sep_sd_400_30$V2, i400s60 = sep_sd_400_60$V2, 
                          i1000j15 = sep_sd_1000_15$V2, i1000j30 = sep_sd_1000_30$V2, i1000j60 = sep_sd_1000_60$V2,
                          i10000j15 = sep_sd_10000_15$V2, i10000j30 = sep_sd_10000_30$V2, i10000j60 = sep_sd_10000_60$V2)
  sep_sd_g3 <- data.frame(i400j15=sep_sd_400_15$V3, i400s30 = sep_sd_400_30$V3, i400s60 = sep_sd_400_60$V3, 
                          i1000j15 = sep_sd_1000_15$V3, i1000j30 = sep_sd_1000_30$V3, i1000j60 = sep_sd_1000_60$V3,
                          i10000j15 = sep_sd_10000_15$V3, i10000j30 = sep_sd_10000_30$V3, i10000j60 = sep_sd_10000_60$V3)
  sep_sd_g4 <- data.frame(i400j15=sep_sd_400_15$V4, i400s30 = sep_sd_400_30$V4, i400s60 = sep_sd_400_60$V4, 
                          i1000j15 = sep_sd_1000_15$V4, i1000j30 = sep_sd_1000_30$V4, i1000j60 = sep_sd_1000_60$V4,
                          i10000j15 = sep_sd_10000_15$V4, i10000j30 = sep_sd_10000_30$V4, i10000j60 = sep_sd_10000_60$V4)
  sep_sd_g5 <- data.frame(i400j15=sep_sd_400_15$V5, i400s30 = sep_sd_400_30$V5, i400s60 = sep_sd_400_60$V5, 
                          i1000j15 = sep_sd_1000_15$V5, i1000j30 = sep_sd_1000_30$V5, i1000j60 = sep_sd_1000_60$V5,
                          i10000j15 = sep_sd_10000_15$V5, i10000j30 = sep_sd_10000_30$V5, i10000j60 = sep_sd_10000_60$V5)
  sep_sd_all <- dplyr::bind_rows(sep_sd_g1,sep_sd_g2,sep_sd_g3,sep_sd_g4,sep_sd_g5)
  #calr
  calr_sd_g1 <- data.frame(i400j15=calr_sd_400_15$V1, i400s30 = calr_sd_400_30$V1, i400s60 = calr_sd_400_60$V1, 
                           i1000j15 = calr_sd_1000_15$V1, i1000j30 = calr_sd_1000_30$V1, i1000j60 = calr_sd_1000_60$V1,
                           i10000j15 = calr_sd_10000_15$V1, i10000j30 = calr_sd_10000_30$V1, i10000j60 = calr_sd_10000_60$V1)
  calr_sd_g2 <- data.frame(i400j15=calr_sd_400_15$V2, i400s30 = calr_sd_400_30$V2, i400s60 = calr_sd_400_60$V2, 
                           i1000j15 = calr_sd_1000_15$V2, i1000j30 = calr_sd_1000_30$V2, i1000j60 = calr_sd_1000_60$V2,
                           i10000j15 = calr_sd_10000_15$V2, i10000j30 = calr_sd_10000_30$V2, i10000j60 = calr_sd_10000_60$V2)
  calr_sd_g3 <- data.frame(i400j15=calr_sd_400_15$V3, i400s30 = calr_sd_400_30$V3, i400s60 = calr_sd_400_60$V3, 
                           i1000j15 = calr_sd_1000_15$V3, i1000j30 = calr_sd_1000_30$V3, i1000j60 = calr_sd_1000_60$V3,
                           i10000j15 = calr_sd_10000_15$V3, i10000j30 = calr_sd_10000_30$V3, i10000j60 = calr_sd_10000_60$V3)
  calr_sd_g4 <- data.frame(i400j15=calr_sd_400_15$V4, i400s30 = calr_sd_400_30$V4, i400s60 = calr_sd_400_60$V4, 
                           i1000j15 = calr_sd_1000_15$V4, i1000j30 = calr_sd_1000_30$V4, i1000j60 = calr_sd_1000_60$V4,
                           i10000j15 = calr_sd_10000_15$V4, i10000j30 = calr_sd_10000_30$V4, i10000j60 = calr_sd_10000_60$V4)
  calr_sd_g5 <- data.frame(i400j15=calr_sd_400_15$V5, i400s30 = calr_sd_400_30$V5, i400s60 = calr_sd_400_60$V5, 
                           i1000j15 = calr_sd_1000_15$V5, i1000j30 = calr_sd_1000_30$V5, i1000j60 = calr_sd_1000_60$V5,
                           i10000j15 = calr_sd_10000_15$V5, i10000j30 = calr_sd_10000_30$V5, i10000j60 = calr_sd_10000_60$V5)
  calr_sd_all <- dplyr::bind_rows(calr_sd_g1,calr_sd_g2,calr_sd_g3,calr_sd_g4,calr_sd_g5)
  #con
  con_sd_g1 <- data.frame(i400j15=con_sd_400_15$V1, i400s30 = con_sd_400_30$V1, i400s60 = con_sd_400_60$V1, 
                          i1000j15 = con_sd_1000_15$V1, i1000j30 = con_sd_1000_30$V1, i1000j60 = con_sd_1000_60$V1,
                          i10000j15 = con_sd_10000_15$V1, i10000j30 = con_sd_10000_30$V1, i10000j60 = con_sd_10000_60$V1)
  con_sd_g2 <- data.frame(i400j15=con_sd_400_15$V2, i400s30 = con_sd_400_30$V2, i400s60 = con_sd_400_60$V2, 
                          i1000j15 = con_sd_1000_15$V2, i1000j30 = con_sd_1000_30$V2, i1000j60 = con_sd_1000_60$V2,
                          i10000j15 = con_sd_10000_15$V2, i10000j30 = con_sd_10000_30$V2, i10000j60 = con_sd_10000_60$V2)
  con_sd_g3 <- data.frame(i400j15=con_sd_400_15$V3, i400s30 = con_sd_400_30$V3, i400s60 = con_sd_400_60$V3, 
                          i1000j15 = con_sd_1000_15$V3, i1000j30 = con_sd_1000_30$V3, i1000j60 = con_sd_1000_60$V3,
                          i10000j15 = con_sd_10000_15$V3, i10000j30 = con_sd_10000_30$V3, i10000j60 = con_sd_10000_60$V3)
  con_sd_g4 <- data.frame(i400j15=con_sd_400_15$V4, i400s30 = con_sd_400_30$V4, i400s60 = con_sd_400_60$V4, 
                          i1000j15 = con_sd_1000_15$V4, i1000j30 = con_sd_1000_30$V4, i1000j60 = con_sd_1000_60$V4,
                          i10000j15 = con_sd_10000_15$V4, i10000j30 = con_sd_10000_30$V4, i10000j60 = con_sd_10000_60$V4)
  con_sd_g5 <- data.frame(i400j15=con_sd_400_15$V5, i400s30 = con_sd_400_30$V5, i400s60 = con_sd_400_60$V5, 
                          i1000j15 = con_sd_1000_15$V5, i1000j30 = con_sd_1000_30$V5, i1000j60 = con_sd_1000_60$V5,
                          i10000j15 = con_sd_10000_15$V5, i10000j30 = con_sd_10000_30$V5, i10000j60 = con_sd_10000_60$V5)
  con_sd_all <- dplyr::bind_rows(con_sd_g1,con_sd_g2,con_sd_g3,con_sd_g4,con_sd_g5)
  # rm
  rm(sep_mean_g1,sep_mean_g2,sep_mean_g3,sep_mean_g4,sep_mean_g5,calr_mean_g1,calr_mean_g2,calr_mean_g3,calr_mean_g4,calr_mean_g5,
     con_mean_g1,con_mean_g2,con_mean_g3,con_mean_g4,con_mean_g5,sep_sd_g1,sep_sd_g2,sep_sd_g3,sep_sd_g4,sep_sd_g5,
     calr_sd_g1,calr_sd_g2,calr_sd_g3,calr_sd_g4,calr_sd_g5,con_sd_g1,con_sd_g2,con_sd_g3,con_sd_g4,con_sd_g5,
     con_mean_400_15,con_mean_400_30,con_mean_400_60,con_mean_1000_15,con_mean_1000_30,con_mean_1000_60,con_mean_10000_15,con_mean_10000_30,con_mean_10000_60,
     calr_mean_400_15,calr_mean_400_30,calr_mean_400_60,calr_mean_1000_15,calr_mean_1000_30,calr_mean_1000_60,calr_mean_10000_15,calr_mean_10000_30,calr_mean_10000_60,
     sep_mean_400_15,sep_mean_400_30,sep_mean_400_60,sep_mean_1000_15,sep_mean_1000_30,sep_mean_1000_60,sep_mean_10000_15,sep_mean_10000_30,sep_mean_10000_60,
     con_sd_400_15,con_sd_400_30,con_sd_400_60,con_sd_1000_15,con_sd_1000_30,con_sd_1000_60,con_sd_10000_15,con_sd_10000_30,con_sd_10000_60,
     calr_sd_400_15,calr_sd_400_30,calr_sd_400_60,calr_sd_1000_15,calr_sd_1000_30,calr_sd_1000_60,calr_sd_10000_15,calr_sd_10000_30,calr_sd_10000_60,
     sep_sd_400_15,sep_sd_400_30,sep_sd_400_60,sep_sd_1000_15,sep_sd_1000_30,sep_sd_1000_60,sep_sd_10000_15,sep_sd_10000_30,sep_sd_10000_60)
  
  cat("tidyr for facet\r")
  # item para
  # すべてをまとめてファセットする
  RMSE_a_res <- data.frame(subject=c(rep(400,300),rep(1000,300),rep(10000,300)) , item=rep(c(rep(15,100),rep(30,100),rep(60,100)),3),
                           SL=tidyr::gather(RMSEa_sep)$value, CC=tidyr::gather(RMSEa_con)$value, calr=tidyr::gather(RMSEa_calr)$value) %>% 
    tidyr::gather(key=method, value=RMSE, -subject, -item)
  
  RMSE_b_res <- data.frame(subject=c(rep(400,300),rep(1000,300),rep(10000,300)) , item=rep(c(rep(15,100),rep(30,100),rep(60,100)),3),
                           SL=tidyr::gather(RMSEb_sep)$value, CC=tidyr::gather(RMSEb_con)$value, calr=tidyr::gather(RMSEb_calr)$value) %>% 
    tidyr::gather(key=method, value=RMSE, -subject, -item)
  # mean
  mean_res <- data.frame(subject=c(rep(400,1500),rep(1000,1500),rep(10000,1500)) , item=rep(c(rep(15,500),rep(30,500),rep(60,500)),3),
                         grade=rep(c(rep("G1",100),rep("G2",100),rep("G3",100),rep("G4",100),rep("G5",100)),3),
                         SL=tidyr::gather(sep_mean_all)$value, 
                         CC=tidyr::gather(con_mean_all)$value, 
                         calr=tidyr::gather(calr_mean_all)$value) %>% 
    tidyr::gather(key=method, value=mean, -subject, -item, -grade)
  # sd
  sd_res <- data.frame(subject=c(rep(400,1500),rep(1000,1500),rep(10000,1500)) , item=rep(c(rep(15,500),rep(30,500),rep(60,500)),3),
                       grade=rep(c(rep("G1",100),rep("G2",100),rep("G3",100),rep("G4",100),rep("G5",100)),3),
                       SL=tidyr::gather(sep_sd_all)$value, 
                       CC=tidyr::gather(con_sd_all)$value, 
                       calr=tidyr::gather(calr_sd_all)$value) %>% 
    tidyr::gather(key=method, value=sd, -subject, -item, -grade)
  # DICC
  dicc_res <- data.frame(condition=condition, item=itemsize, subject=subjectsize,calr=calr_dicc, SL=sep_dicc, con=con_dicc) %>% 
    tidyr::gather(key=method, value=DICC, -item, -subject, -condition)
  # DICC_pdist
  dicc_pdist_res <- data.frame(condition=condition, item=itemsize, subject=subjectsize,calr=calr_dicc_pdist, SL=sep_dicc_pdist, con=con_dicc_pdist) %>% 
    tidyr::gather(key=method, value=DICC, -item, -subject, -condition)
  
  # まずは全実験結果を一覧でプロット
  
  leg_lab <- c("A1","A2","A3","B1","B2","B3","C1","C2","C3") # legend
  leg <- c("A1: 受検者数10,000　項目数15","A2: 受検者数10,000　項目数30","A3: 受検者数10,000　項目数60",
           "B1: 受検者数1,000　項目数15","B2: 受検者数1,000　項目数30","B3: 受検者数1,000　項目数60",
           "C1: 受検者数400　項目数15","C2: 受検者数400　項目数30","C3: 受検者数400　項目数60") # legend
  
  calrplot_a <- RMSEa_calr %>% tidyr::gather(key = condition, value = RMSE) %>% 
    ggplot(aes(x=condition,y=RMSE,fill=condition))+
    ylim(0,0.3) + 
    #geom_boxplot()+ 
    geom_violin()+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) +
    scale_fill_manual(values = c("#1874CD","#1874CD","#1874CD","#EE6363","#EE6363","#EE6363","#008B45","#008B45","#008B45"),
                      labels = leg) + 
    scale_x_discrete(labels = leg_lab)
  
  calrplot_b <- RMSEb_calr %>% tidyr::gather(key = condition, value = RMSE) %>% 
    ggplot(aes(x=condition,y=RMSE,fill=condition))+
    #ylim(0,2) + 
    #geom_boxplot()+ 
    geom_violin()+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) +
    scale_fill_manual(values = c("#1874CD","#1874CD","#1874CD","#EE6363","#EE6363","#EE6363","#008B45","#008B45","#008B45"),
                      labels = leg) + 
    scale_x_discrete(labels = leg_lab)
  
  sepplot_a <- RMSEa_sep %>% tidyr::gather(key = condition, value = RMSE) %>% 
    ggplot(aes(x=condition,y=RMSE,fill=condition))+
    ylim(0,0.3) + 
    #geom_boxplot()+ 
    geom_violin()+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) +
    scale_fill_manual(values = c("#1874CD","#1874CD","#1874CD","#EE6363","#EE6363","#EE6363","#008B45","#008B45","#008B45"),
                      labels = leg) + 
    scale_x_discrete(labels = leg_lab)
  
  sepplot_b <- RMSEb_sep %>% tidyr::gather(key = condition, value = RMSE) %>% 
    ggplot(aes(x=condition,y=RMSE,fill=condition))+
    #ylim(0,1) + 
    #geom_boxplot()+ 
    geom_violin()+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) +
    scale_fill_manual(values = c("#1874CD","#1874CD","#1874CD","#EE6363","#EE6363","#EE6363","#008B45","#008B45","#008B45"),
                      labels = leg) + 
    scale_x_discrete(labels = leg_lab)
  
  conplot_a <- RMSEa_con %>% tidyr::gather(key = condition, value = RMSE) %>% 
    ggplot(aes(x=condition,y=RMSE,fill=condition))+
    #ylim(0,0.3) + 
    #geom_boxplot()+ 
    geom_violin()+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) +
    scale_fill_manual(values = c("#1874CD","#1874CD","#1874CD","#EE6363","#EE6363","#EE6363","#008B45","#008B45","#008B45"),
                      labels = leg) + 
    scale_x_discrete(labels = leg_lab)
  
  conplot_b <- RMSEb_con %>% tidyr::gather(key = condition, value = RMSE) %>% 
    ggplot(aes(x=condition,y=RMSE,fill=condition))+
    #ylim(0,2) + 
    #geom_boxplot()+ 
    geom_violin()+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=15,face="bold")) +
    scale_fill_manual(values = c("#1874CD","#1874CD","#1874CD","#EE6363","#EE6363","#EE6363","#008B45","#008B45","#008B45"),
                      labels = leg) + 
    scale_x_discrete(labels = leg_lab)
  
  # ggplot_facet
  RMSE_a_plot <- ggplot(RMSE_a_res, aes(y=RMSE, x=method, fill=method, colour=method)) + 
    #geom_boxplot()+ 
    geom_violin()+
    #ylim(0,0.4) +
    facet_grid(item~subject) +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
    scale_fill_manual(values = c(cols[1],cols[5],cols[3])) +# 色分けを独自定義
    scale_colour_manual(values = c(cols[2],cols[6],cols[4])) # 色分けを独自定義
  
  RMSE_b_plot <- ggplot(RMSE_b_res, aes(y=RMSE, x=method, fill=method, colour=method)) + 
    #geom_boxplot()+ 
    geom_violin()+
    #ylim(0,2) + # 外れ値がいくつか消される。
    facet_grid(item~subject) +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
    scale_fill_manual(values = c(cols[1],cols[5],cols[3])) +# 色分けを独自定義
    scale_colour_manual(values = c(cols[2],cols[6],cols[4])) # 色分けを独自定義
  
  # DICC
  # DICC_nomal
  dicc_plot <- ggplot(dicc_res, aes(y=DICC, x=method, fill=method)) + 
    geom_bar(stat = "identity") + 
    facet_grid(item~subject) +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
    scale_fill_manual(values = c(cols[2],cols[6],cols[4]), lab=c("calr","CC","SL"))
  #scale_fill_manual(values = c("#1874CD","#EE6363","#008B45")) # 色分けを独自定義(カラフル)
  #scale_fill_manual(values = c("#030303","#454545","#ADADAD")) # モノクロ（論文，白黒用）
  
  # DICC_pdist
  dicc_pdist_plot <- ggplot(dicc_pdist_res, aes(y=DICC, x=method, fill=method)) + 
    geom_bar(stat = "identity") + 
    facet_grid(item~subject) +
    theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
    scale_fill_manual(values = c(cols[2],cols[6],cols[4]), lab=c("calr","CC","SL")) # 色分けを独自定義(カラフル)
  
  mean_res_plot <- mean_res %>% ggplot(aes(x=grade,y=mean,fill=method,colour=method)) + 
    geom_boxplot()+ 
    #geom_violin()+
    #ylim(0,2)+
    facet_grid(item~subject) + 
    #geom_hline(yintercept = c(0,0.4,0.8,1.2,1.6), colour="blue", linetype="dashed") + 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
    scale_fill_manual(values = c(cols[1],cols[5],cols[3]), lab=c("calr","CC","SL")) +# 色分けを独自定義
    scale_colour_manual(values = c(cols[2],cols[6],cols[4])) # 色分けを独自定義
  
  sd_res_plot <- sd_res %>% ggplot(aes(x=grade,y=sd,fill=method,colour=method)) + 
    geom_boxplot()+ 
    #geom_violin()+
    facet_grid(item~subject) + 
    theme(axis.text=element_text(size=12),axis.title=element_text(size=12,face="bold")) +
    scale_fill_manual(values = c(cols[1],cols[5],cols[3]), lab=c("calr","CC","SL")) +# 色分けを独自定義
    scale_colour_manual(values = c(cols[2],cols[6],cols[4])) # 色分けを独自定義
  #geom_hline(yintercept = c(1,0.9,0.8,0.7,0.6), colour="blue", linetype="dashed")
  

  
  cat("output plot as meta file\r")
  
  setwd(paste0("/Users/sep10_z1vk0al/OneDrive/Documents/00_修士論文/plot-for-M/4-1_plot_assy/",sim))
  
  # グラフ出力の前に，記述統計量を確認
  sink(file="tapply_stat.txt")
  cat("RMSE_a\n")
  with(RMSE_a_res,tapply(RMSE, list(method,subject,item), mean, simplify=F))
  cat("RMSE_b\n")
  with(RMSE_b_res,tapply(RMSE, list(method,subject,item), mean))
  cat("pdist_mean\n")
  with(mean_res, tapply(mean, list(method,grade,subject,item), mean))
  cat("pdist_sd\n")
  with(sd_res, tapply(sd, list(method,grade,subject,item), mean))
  sink()
  
  emf(paste0(sim,"_sepplot_a.emf"), width = 9, height = 5)
  sepplot_a
  dev.off()
  
  emf(paste0(sim,"_sepplot_b.emf"), width = 9, height = 5)
  sepplot_b
  dev.off()
  
  emf(paste0(sim,"_calrplot_a.emf"), width = 9, height = 5)
  calrplot_a
  dev.off()
  
  emf(paste0(sim,"_calrplot_b.emf"), width = 9, height = 5)
  calrplot_b
  dev.off()
  
  emf(paste0(sim,"_conplot_a.emf"), width = 9, height = 5)
  conplot_a
  dev.off()
  
  emf(paste0(sim,"_conplot_b.emf"), width = 9, height = 5)
  conplot_b
  dev.off()
  
  emf(paste0(sim,"_RMSE_a_facet.emf"), width=7, height = 5)
  RMSE_a_plot
  dev.off()
  
  emf(paste0(sim,"_RMSE_a_facet_zoom.emf"), width=7, height = 5)
  RMSE_a_plot+ylim(0,0.25)
  dev.off()
  
  emf(paste0(sim,"_RMSE_b_facet.emf"), width=7, height = 5)
  RMSE_b_plot
  dev.off()
  
  emf(paste0(sim,"_RMSE_b_facet_zoom.emf"), width=7, height = 5)
  RMSE_b_plot+ylim(0,0.5)
  dev.off()
  
  emf(paste0(sim,"_DICC_facet_plot.emf"), width = 7,height = 5)
  dicc_plot
  dev.off()
  
  emf(paste0(sim,"_DICC_pdist_plot.emf"), width=7, height=5)
  dicc_pdist_plot
  dev.off()
  
  emf(paste0(sim,"_mean_plot.emf"), width=8, height=5)
  mean_res_plot+
    geom_hline(yintercept = 0, linetype="dotdash")+
    geom_hline(yintercept = 0.4, linetype="dotdash")+
    geom_hline(yintercept = 0.8, linetype="dotdash")+
    geom_hline(yintercept = -0.4, linetype="dotdash")+
    geom_hline(yintercept = -0.8, linetype="dotdash")
  
  dev.off()
  
  emf(paste0(sim,"_sd_plot.emf"), width=8, height=5)
  sd_res_plot+
    geom_hline(yintercept = 1, linetype="dotdash")
  dev.off()
  
} # end of sim

# 以下，シミュレーションの説明のための作図

emf("simplot.emf", width = 7, height = 5)
curve(1/(1+exp(-1.7*1*(x-0))), -4, 4, ylab = "正答確率", xlab = "受検者能力値" , sub = "項目パラメタ:a=1, b=0")
abline(v=1, lty=2)
arrows(1,0,1,1/(1+exp(-1.7*1*(1-0))), col = 2, length = 0.1, code = 3)
arrows(1,1,1,1/(1+exp(-1.7*1*(1-0))), col = 4, length = 0.1, code = 3)
dev.off()

#the graph of item characteristic curve
twoplm <- function(theta,a,b,c){
  c+(1-c)/(1+exp(-1.702*a*(theta-b)))
}

theta <- data.frame(theta=seq(seq(-6, 6, 100)))
theta <- ggplot(theta, aes(x=theta))
zz <- theta + stat_function(fun=twoplm, args=list(a=1,b=0,c=0))+ 
  labs(title="2PLMの項目反応関数",x="受検者能力値",y="正答確率")+
  scale_x_continuous(breaks=seq(-4,4,by=1),limits=c(-4,4)) + 
  geom_vline(xintercept = 1, lty = 2)
zz


# 能力分布

#theta <- data.frame(r=rnorm(10000,0,1))
theta <- data.frame(theta=seq(seq(-4, 4, 100)))

g <- ggplot(theta, aes(x=theta))
#geom_density()
for(j in 1:5){
  g <- g + stat_function(fun=dnorm, args = list(mean=thetaMS[j,1], sd=thetaMS[j,2]),
                         colour = j+1, size = 1)
  # geom="area", fill = "color"で塗りつぶしは可能
} 
g <- g +  scale_x_continuous(breaks = seq(-4,4, by=1), limits = c(-4,4)) +
  labs(x=TeX("$\\theta$"), y=TeX("$P(\\theta)$"))
g

emf("thetadist.emf", width=7, height = 5)
dev.off()


gg <- ggplot(theta, aes(x=theta)) + 
  stat_function(fun=dlnorm, args = list(mean=0, sd=0.25), geom="area", fill="#EE6363") + 
  labs(x="a",y="Probability density") +
  scale_x_continuous(limits = c(0,2))
emf("discrimidist.emf", width=7, height = 5)
gg
dev.off()

ggg <- ggplot(theta, aes(x=theta)) + 
  stat_function(fun=dnorm, args = list(mean=0.8, sd=1.5), geom="area", fill="#1874CD") + 
  labs(x="b(st)",y="Probability density") +
  scale_x_continuous(limits = c(-5.2,6.8))
emf("diffdist.emf", width=7, height = 5)
ggg
dev.off()

