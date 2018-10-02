# myfanction and useful function


#---------------------------------------------------------------------------#
#MML-EM関数
#---------------------------------------------------------------------------#

# Rcpp version
library(Rcpp)

sourceCpp("C:/Users/sep10_z1vk0al/OneDrive/Documents/12_R/MML-EM/MML-EM/MML-EM_debug.cpp")
sourceCpp("~/OneDrive/Documents/12_R/MML-EM/MML-EM/MML-EM_debug.cpp")


MML_EM <- function(x, gh=TRUE, e1=0.001, e2=0.0001, D=1.702, method="FS", gr=dtj
                   , iter=200, maxinitial = 20, fc=3, lc=ncol(x)
                   , N = 31, maxtheta = 4.75, mintheta = -4.75
                   , Fname="data", output = T, FA = F){
  
  #-------------------------------------------------------------------------#
  #
  # method  :Select which one. "NR" is Newton-Raphson, "FS" is Fisher-Scoring,
  #          "BFGS" is Qusai-Newton.
  #
  # output csv file is adaptive for EasyEstimation, Ryuichi Kumagai (2009).
  #
  # 
  # この関数はRのパッケージである，dplyr, pracma, psychに依存しています。関数を走らせる前に，これらのパッケージをインストール&ライブラリしておいてください。
  # 
  #
  # 引数について
  #
  # x		: 項目反応データ。デフォルトの設定では1行目にID，2行目に母集団変数が入っていることを想定しています。
  # gh		: 区分求積法にガウス・エルミート法を使うかどうか。Fにすると，等区分求積法で実行します。
  # e1		: EMステップの収束基準です。項目パラメタの変化量で決めています。
  # e2		: Mステップの収束基準です。同じく，項目パラメタの変化量で決めています。
  # D		: 尺度因子，定数です。
  # method	: Mステップの最大化の方法です。デフォルトはFS＝フィッシャースコアリングですが，ほかにもNR＝ニュートンラフソン法や，Rのoptim関数の諸手法を利用できます。
  # 		: FSは初期値に敏感な手法のようです。NRは，逆行列の計算でつまづくことがあります。BFGSやCGは収束基準に向かって漸近的に最大化しないことがあります。
  # gr		: グラディエント指定用の関数です。使いません。
  # iter	: EMステップの最大繰り返し数です。デフォルトは200回です。
  # maxinitial	: 初期値の最大値についてです。この数値以上の初期値が検出された場合，識別力の初期値を1，困難度の初期値を0にして，計算を実行します。
  # fc		: 反応データの開始列です。デフォルトは3です。
  # lc		: 反応データの終了列です。デフォルトは最終列です。
  # N		: 区分求積の分点数です。デフォルトは31です。
  # maxtheta	: 区分求積の最大値です。デフォルトは;4.75です。
  # mintheta	: 区分求積の最小値です。デフォルトは-4.75です。
  # Fname	: 出力ファイル名です。
  # output	: csvファイルを書き出すかどうかです。Fにすれば出力しません。
  # 
  #
  #
  #-------------------------------------------------------------------------#
  
  library(psych, dplyr, pracma)
  message("Reading item response data; 項目反応データを読み込んでいます。")
  xall <- x[,fc:lc]
  
  ni <- nrow(xall)
  message("the number of subjects is ",ni)
  nj <- ncol(xall)
  message("the number of items is ",nj)
  
  
  status <- c("NULL")
  t0 <- matrix(0,nj,2)
  initial <- matrix(0,nj,2)
  r <- numeric(0)
  p <- numeric(0)
  EVP <- numeric(0)
  
  if(FA == T){
    message(paste("Checking reduction of eigen values based on tetrachoric correlation coefficient;テトラコリック相関係数にもとづく固有値の減衰状況を確認します。"))
    
    
    
    CHECK <- try({tetrachoric(xall)})
    
    if(class(CHECK) != "try-error"){
      FAdata <- fa(xall,nfactors = 1,rotate = "none",cor="tet")
    }else{
      message("\nテトラコリック相関係数の計算に失敗しました")
    }
    
    if(class(CHECK) == "try-error"){
      message("点双列相関係数と通過率にもとづいて初期値を計算します。")
      status <- c("Failuer FA")
      r <- as.vector(cor(rowSums(xall,na.rm = T),xall,use="pairwise.complete.obs"))
      p <- colMeans(xall,na.rm = T)
      t0 <- matrix(nrow=nj,ncol=2)
      t0[,1] <- D*r/sqrt(1-r^2)
      t0[,2] <- -log(p/(1-p))  #log odds ratio
    } else if(min(FAdata$loadings)<0){
      message("\n1因子での因子分析の結果，負の因子負荷量が得られました。")
      message("点双列相関係数と通過率にもとづいて初期値を計算します。")
      status <- c("FA")
      r <- as.vector(cor(rowSums(xall,na.rm = T),xall,use="pairwise.complete.obs"))
      p <- colMeans(xall,na.rm = T)
      t0 <- matrix(nrow=nj,ncol=2)
      t0[,1] <- D*r/sqrt(1-r^2)
      t0[,2] <- -log(p/(1-p))  #log odds ratio
    }else{  
      message(paste(round(FAdata$values[1],4),round(FAdata$values[2],4),round(FAdata$values[3],4)),"\n")
      EVP <- plot(FAdata$values, lty=1,type="o",xlab="compnent number", ylab="eigen values", main="scree plot")
      status <- c("FA")
      message("\nテトラコリック相関係数と通過率にもとづいて初期値を計算します。")
      r <- FAdata$loadings
      p <- colMeans(xall,na.rm = T)
      t0 <- matrix(nrow=nj,ncol=2)
      t0[,1] <- D*r/sqrt(1-r^2)
      t0[,2] <- qnorm(p,0,1,lower.tail=FALSE)/r # こっちの初期値でもいける。
    }
  } else {
    status <- c("Failuer FA")
    r <- as.vector(cor(rowSums(xall,na.rm = T),xall,use="pairwise.complete.obs"))
    p <- colMeans(xall,na.rm = T)
    t0 <- matrix(nrow=nj,ncol=2)
    t0[,1] <- D*r/sqrt(1-r^2)
    t0[,2] <- -log(p/(1-p))  #log odds ratio
  }



  
  
  
  message("初期値の最大絶対値",round(max(abs(t0)),4))
  if(max(abs(t0)) > maxinitial){
    message("初期値の絶対値が",maxinitial,"以上でした。不適切な初期値の可能性があります。")
    message("不適切な初期値を識別力1，困難度0に設定して計算します。")
    t0[abs(t0[,1])>maxinitial,1] <- 1
    t0[abs(t0[,2])>maxinitial,2] <- 0
  }
  
  initial <- t0
  
  #Check the initial value
  #message("a-parameter")
  #for (j in 1:nj){
  #  item <- t0[j,1]
  #  message(round(item, digits = 3))
  #}
  #message("b-parameter")
  #for (j in 1:nj){
  #  item <- t0[j,2]
  #  message(round(item, digits = 3))
  #}
  #message(lapply(t0,round, digits=3))
  
  if(gh == TRUE){   #ghにガウスエルミートのファイルが指定されていた場合
    message("ガウスエルミート求積法が指定されました。")
    message("分点数は",N ,"です。")
    nm <- N
    Xm <- gaussHermite(N)$x*sqrt(2)
    Wm <- gaussHermite(N)$w/sqrt(pi)  #caribration
  }else{
    message("等区分求積法が指定されました。")
    message("分点数は",N ,"です。最大の特性値は",maxtheta,"で，最小の特性値は",mintheta,"です。")
    nm <- N
    Xm <- seq(mintheta, maxtheta, length.out = N)
    Wm <- dnorm(Xm, 0, 1)
    Wm <- Wm/sum(Wm)
  }
  
  ptheta <- function(theta, a, b, D){
    1/(1+exp(-D*a*(theta-b)))
  }
  L <- function(xi,theta, a, b, D){
    exp(sum(log(ptheta(theta,a,b,D))*xi+log(1-ptheta(theta,a,b,D))*(1-xi),na.rm = T))
  }
  ELL <- function(theta,param,r,N,D){
    a <- param[1]
    b <- param[2]
    sum(r*log(ptheta(theta,a,b,D))+(N-r)*log(1-ptheta(theta,a,b,D)),na.rm = T)
  }
  dtj <- function(theta,param,r,N,D){
    a <- param[1]
    b <- param[2]
    da <- D*sum((theta-b)*(r-N*ptheta(theta,a,b,D)),na.rm = T)
    db <- -D*a*sum(r-N*ptheta(theta,a,b,D),na.rm = T)
    c(da,db)
  }
  
  #Hessian matrix for Nweton-Rapthon Algorithm
  Htj <- function(theta,param,r,N,D){
    a <- param[1]
    b <- param[2]
    d2aa <- -D^2*sum((theta - b)^2*N*ptheta(theta,a,b,D)*(1-ptheta(theta,a,b,D)),na.rm = T)
    d2bb <- -D^2*a^2*sum(N*ptheta(theta,a,b,D)*(1-ptheta(theta,a,b,D)),na.rm = T)
    d2ab <- -D*sum(r-N*ptheta(theta,a,b,D),na.rm = T)
    matrix(c(d2aa,d2ab,d2ab,d2bb),nrow=2)
  }
  
  #Fisher's Information matrix of item j
  Itj <- function(theta,param,N,D){
    a <- param[1]
    b <- param[2]
    daa <- sum(N*ptheta(theta,a,b,D)*(1-ptheta(theta,a,b,D))*(theta-b)^2,na.rm = T)
    dab <- sum(-N*ptheta(theta,a,b,D)*(1-ptheta(theta,a,b,D))*(a*(theta-b)),na.rm = T)
    dbb <- sum(N*ptheta(theta,a,b,D)*(1-ptheta(theta,a,b,D))*a^2,na.rm = T)
    (D^2)*matrix(c(daa,dab,dab,dbb),nrow=2)
  }
  
  #convergence criteria
  conv1 <- 0;conv2 <- 0
  
  #n of iteraiton 
  count1 <- 0
  
  message("数値最適化の手法は",method,"が選択されました。")
  
  if(method == "SANN" || method == "Nelder-Mead") gr <-  NULL
  
  system.time(
    while(conv1==0){
      count1 <- count1+1
      #E step starts here
      message(paste("第",count1,"回目のEステップを開始します。"))
      
      #estimate likelihood: L(parameter, nodes | response)
      Flim_apply <-function(xall,a,b,Xm,D){
        L <- function(xi,theta, a, b, D){
          exp(sum(log(ptheta(theta,a,b,D))*xi+log(1-ptheta(theta,a,b,D))*(1-xi),na.rm = T))
        }
        Lim_i <- apply(xall,1,L,a=a,b=b,theta=Xm,D=D)
      }
      #preparation empty matrix
      Lim_apply <- matrix(nrow = ni, ncol = nm)
      #Lim is a matrix, the likelihood of subject i for item j in [i,j] element
      Lim_apply <- apply(matrix(Xm, ncol = 1),1,Flim_apply,xall=xall,a=t0[,1],b=t0[,2],D=D)
      
      #calculate the weight for a posteriori distribution
      #product likelihood and weight of a priori distribution
      Fgim_apply <- function(Lim_i,Wm){
        Lim_i*Wm/sum(Lim_i*Wm)
      }
      Gim_apply <- apply(Lim_apply,1,Fgim_apply,Wm=Wm)
      Gim_apply <- t(Gim_apply)
      
      #Expectation of subject for each nodes
      #Nm <- colSums(Gim_apply,na.rm = T)
      
      Nm <- matrix(ncol = N, nrow = nj)

    	for(i in 1:nj){　#　欠測値を考慮した，期待度数の計算
	      resp <- xall[,i]
	      sub <- !is.na(resp)

    	  Fnm <- function(Gim_j, sub){
	        sum(Gim_j * sub, na.rm = T)
	      }

	      Nm[i,] <- apply(Gim_apply, 2, Fnm, sub = sub)
	    }

      xall[is.na(xall)] <- 0　#　欠測値を考慮した正答期待度数の計算
      rjm <- t(xall)%*%Gim_apply
      xall <- x[,fc:lc]
      
      #M step starts here
      message("第",count1,"回目のMステップを開始します。")
      t0m <- t1m <- t0
      
      if(method == "NR"){
        for (j in 1:nj){
          conv2<- 0
          while(conv2 == 0){
            (t1m[j,] <- t0m[j,] - solve(Htj(Xm,t0m[j,],rjm[j,],Nm[j,],D))%*%dtj(Xm,t0m[j,],rjm[j,],Nm[j,],D))
            if (max(abs(t1m[j,]-t0m[j,])) < e2) conv2 <- 1
            t0m[j,] <- t1m[j,]
          }
        }
      }else if(method == "FS"){
        for (j in 1:nj){
          conv2<- 0
          while(conv2==0){
            t1m[j,] <- t0m[j,] + solve(Itj(Xm,t0m[j,],Nm[j,],D))%*%dtj(Xm,t0m[j,],rjm[j,],Nm[j,],D)
            if (all(abs(t1m[j,]-t0m[j,]) < e2)) conv2 <- 1
            t0m[j,] <- t1m[j,]
          }
        }
      }else{
        for(j in 1:nj){
          BFGS <- optim(par=t0[j,], fn=ELL,theta=Xm, r=rjm[j,], N=Nm[j,], D=D,method=method, gr=gr,
                        control = list(reltol = e2,fnscale = -1,trace=F))#reltol = e2,#gr=dtj,method="BFGS"
          t1m[j,] <- BFGS$par
          #------------------------------------------------------------------------------------------------#
          #R function for Optimization
          #"BFGS"だと収束しなかった。CGは180秒近くかかって収束が遅い。"L-BFGS-B"だとCGよりも少し早く収束した。
          #デフォルトの方"Nelder-Mead"法もなかなか早い。80秒を切っている。SANNは遅すぎて使い物にならなそう。
          #------------------------------------------------------------------------------------------------#
        }
        
      }
      
      if((max(abs(t1m-t0))<e1)||(count1==iter)){
        conv1 <- 1
        message("項目パラメタの計算が終了しました。")
        para <- round(t1m,5)
        message("推定の標準誤差を計算します。")
        Lim_apply <- apply(matrix(Xm, ncol = 1),1,Flim_apply,xall=xall,a=t1m[,1],b=t1m[,2],D=D)
        Gim_apply <- apply(Lim_apply,1,Fgim_apply,Wm=Wm)
        Gim_apply <- t(Gim_apply)
        Nm <- matrix(ncol = N, nrow = nj)

	for(i in 1:nj){
  
	  resp <- xall[,i]
	  sub <- !is.na(resp)
  
	  Fnm <- function(Gim_j, sub){
    
	    sum(Gim_j * sub, na.rm = T)
    
	  }
  
	  Nm[i,] <- apply(Gim_apply, 2, Fnm, sub = sub)
  
  
	}
        SE <- matrix(0,nrow=nj,ncol=2);colnames(SE) <- c("a","b")
        #      itj <- apply(t1m,1,Itj,theta=Nm,N=Xm)
        try(
          for(j in 1:nj){
            SE[j,] <- round(sqrt(diag(solve(Itj(theta=Xm,t1m[j,],N=Nm[j,],D)))),5)
          }
        )
        message("推定の標準誤差の計算が終了しました。")
        #推定母集団分布
        prob <- as.matrix(Nm/sum(Nm), ncol = 1)
        theta <- as.matrix(Xm, ncol=1)
        population_mean <- mean(prob%*%theta)
        pmat <- matrix(population_mean,nrow(prob),ncol(prob),byrow=T)
        thetamat <- matrix(theta,nrow(prob),ncol(prob),byrow=F)
        population_vari <- mean(colSums(((thetamat - pmat)^2)*prob))
        EPD <- data.frame(Mean = population_mean, SD = sqrt(population_vari))
        message(("推定母集団分布の計算が終了しました。"))
        #項目パラメタファイルの体裁を整える
        ItemID <-  as.matrix(paste0(" Item",formatC(c(1:nj),flag="0",width="3")))
        PL <- matrix(" 2PL",nj,1)
        c <- matrix(0,nj,1)
        param <- cbind(ItemID,PL,round(t1m,5),c)
        param <- as.data.frame(param);colnames(param) <- c("ItemID","model","slope","location","asymptote")
        if(status != "FA"){FAdata <- list (values=numeric(0), loadings=numeric(0))}
        item_para_result <- list(Method=method ,D=paste0("D = ",D),param=param, StandardError = round(SE,5), EstimatedPopulationDistribution = data.frame(theta=theta,prob=colMeans(prob))
                                 , EstimatedPopulationDistribution = EPD, iteration=count1
                                 , InitialValue=round(data.frame(a=initial[,1],b=initial[,2]),5),ScreePlot=EVP, EigenValues=FAdata$values, Lordings=FAdata$loadings
        )

        if(output == T){
          
          message("リザルトファイルを書き出しています。")
          
          try({
            for (n in 1:length(item_para_result)){
              if(n == 1) {
                write.table(item_para_result[1],file=paste0(Fname,"_ParaResult.csv"),sep = ",",quote = F, row.names = F)
                write.table(item_para_result[2],file=paste0(Fname,"_Para.csv"),sep = "",quote = F, col.names = F, row.names = F)
                write.table(item_para_result[3],file=paste0(Fname,"_Para.csv"),sep = ",",quote = F, append = T, col.names = F, row.names = F)
              }else if(n >1){
                write.table(item_para_result[n],file=paste0(Fname,"_ParaResult.csv"),sep = ",",quote = F, append = T, row.names = F)
                write.table("\n",file=paste0(Fname,"_ParaResult.csv"),sep = "",quote = F, append = T, col.names = F, row.names = F)
              }
            }
          })
        }
        
        return(item_para_result)
        
      } else{
        message("EMアルゴリズムは第",count1,"回の反復計算が終了しました。")
        message("周辺対数尤度は",sum(log(Lim_apply%*%Wm)),"です。")
        message("項目パラメタの最大変化量は",max(abs(t1m-t0)),"です。\n")
        t0 <- t1m
      }
    }
  )
}




#---------------------------------------------------------------------------#
#能力パラメタ，推算値推定関数
#---------------------------------------------------------------------------#


estheta <- function(xall, param, est, nofrands=10, method="NR", file="default", output=FALSE, IDc=1, Gc=2, ITEMc=3,
                    gh = TRUE, N = 31, maxtheta = 6, mintheta = -6, mu=0, sigma=1){
  
  
  #--------------------------------------------#
  # This function was edited 
  # 		by D. EJIRI and T. SHIBUYA.
  # 		at Tohoku University in Japan.
  # 		(DATE)
  # 		The Original syntax was written by T. SHIBAYAMA.
  #
  # Oct. 10 2017    元シンタックスの修正（江尻）
  # Oct. 22 2017    デバッグ（澁谷）
  # Nov. 10 2017    結果ファイルをCSVデータとして出力するように設定（江尻）
  # Feb. 21 2018    関数内でMAP推定をおこない，その推定値を事後分布の最大確率密度に設定するように変更（澁谷）
  #                 推定アルゴリズムは基本的にフィッシャースコアリングを用い，100回以上の反復計算をおこなっても収束しない場合には二分法を用いるようにした。
  # Apr.  3 2018    計算の高速化のために，forによるループを回避してapply関数を用いるように変更した（澁谷）
  # Apr.  5 2018    集団統計量の補間の標準誤差を計算する関数を追加（澁谷）
  # May. 18 2018    推算値推定関数を拡張し，MLE，MAP，EAPを推定可能な関数として更新。関数名をplausibleからesthetaに変更。
  #                 WLE推定法は計算結果の検算ができていないが，おおよそ妥当そうな値が計算できている。　※まだ実用段階ではない。
  # Sep. 10 2018    3PLに拡張，一部の反応パタンにおいて最尤推定が失敗するバグを修正。
  #
  #
  #------------------------------引数の仕様について-----------------------------#
  #
  # xall     : 2値の項目反応データ。デフォルトではIDが一列目，母集団変数が2列目，項目反応パタンが３列目から最終行まで入力されているデータを想定。
  #          : 母集団変数が不要な場合はGc=0と入力する。
  # para     : 項目反応データ.一列目には識別力を，二列目には困難度が入っているデータ。
  # nofrands : 発生させたい推算値の組数
  # gh       : ガウスエルミート求積法に関するオプション。FALSEならば等区分求積法で計算は実行される。
  #          : 分点と重みの計算にはpracmaパッケージのgaussHermite関数を利用した。
  # N        : 区分求積の分点数。
  # mintheta : 等区分求積法のθの最小値。
  # maxtheta : 等区分求積法のθの最大値。
  # output   : 推算値の計算結果をCSVとして出力するかどうか。デフォルトはFALSE。
  # file     : CSVファイル名。ただし拡張子は除く。デフォルトは"default"
  #
  #  ～推算値の計算手順～
  # 1）EAP推定関数の定義，EAP推定値と周辺確率の計算
  # 2）MAP推定関数の定義，MAP推定値の計算
  # 3）フォンノイマンの棄却法にもとづいて，推算値を発生させる。
  #
  #--------------------------------------------#
  
  message("データチェック")
  ID <- xall[,IDc]
  if(Gc == 0){
    group <- rep(1, nrow(xall))
    G <- 1
    x.all <- xall[,ITEMc:ncol(xall)]
  }else{
    group <- xall[,Gc]
    G <- max(as.numeric(group))
    if(is.vector(xall)){
      x.all <- xall[ITEMc:length(xall)]
      x.all <- matrix(x.all,nrow=1)
    }
    x.all <- xall[,ITEMc:ncol(xall)]
  }
  #remove a=0 item parameter and response column
  message("識別力が0の項目を削除します。")
  a <- param[,1]
  x.all <- x.all[, a != 0]
  message(sum((a == 0)*1),"個の項目が削除されました")
  
  #  message("全回答がNAの受験者の項目反応データを削除します。")
  #  x.all <- x.all[!is.na(apply(x.all, 1, sum)), ]
  # set item parameter data
  param <- param[param[,1] != 0, ]
  a <- param[,1]
  b <- param[,2]
  c <- param[,3]
  
  # Number of Subjects"
  n <- nrow(x.all)
  ng <- numeric(G)
  # number of items
  m <-length(a)
  xscore <- rowSums(x.all,na.rm=T)
  message("IDと項目反応データ，項目パラメタを確認しました。")
  message("母集団",G)
  message("受検者数",n,"人")
  message("項目数",m,"個でした")
  message(est,"の推定を開始します。")
  
  # 集団ごとのテスト項目数を確認
  groupitem <- numeric(G)
  for(g in 1:G){
    key <- group==g
    temp <- x.all[key,]
    temp <- temp %>% colSums(na.rm = T) 
    groupitem[g] <- temp[temp!=0] %>% length()
  }
  #cat(groupitem)
  
  if((est == "PVs") || (est == "EAP")){
    if(gh == TRUE){   #ガウスエルミートで区分求積を実行する場合
      npoint <- N
      message("ガウスエルミート求積法が指定されました。")
      message("分点数は",npoint,"です。")
      gh <- pracma::gaussHermite(npoint)
      xnodes <- gh$x*sqrt(2)
      weight <- gh$w/sqrt(pi)  #caribration
      ###  Nodes and weights for the npoint Gauss-Hermite formula ####
      #
      #   Integrate of f(x) * exp(-x^2) = Sum of f(x) * weight
      #
      #      ex. Integrate of exp(-x^2) = sqrt(pi) = 1.772454
      #          Sum of weights = 1.772454
      #--------------------------------------------------------------#
    }else{
      message("等区分求積法が指定されました。")
      message("分点数は",N,"です。最大の特性値は",maxtheta,"で，最小の特性値は",mintheta,"です。")
      npoint <- N
      xnodes <- seq(mintheta, maxtheta, length.out = N)
      weight <- dnorm(xnodes, mu, sigma)
      weight <- weight/sum(xnodes)
    } 
  }
  
  #----------------------------------------------------------#
  #EAP estimater
  #----------------------------------------------------------#
  
  
  FthetaEAP <- function(xi,theta,w,a,b,c){
    #--------------------------------------------#
    #MAP推定をおこなう関数。
    #xi    : 受検者iの項目反応データ
    #theta : 事前分布の分点
    #w     : 事前分布の分点に対応する重み
    #a&b   : パラメタファイル。aが識別力，bが困難度
    #--------------------------------------------#
    theta <- as.matrix(theta,byrow=T)
    xi <- as.numeric(xi)
    
    # P(theta) in two-parameter logisticmodel
    ptheta <- function(theta,a,b,c){
      D <- 1.702
      c+(1-c)/(1+exp(-D*a*(theta-b)))
    }
    
    # 対数尤度
    LL <- function(u,theta,a,b,c){
      sum(u*log(ptheta(theta,a,b,c))+(1-u)*log(1-ptheta(theta,a,b,c)),na.rm = T)
    }
    
    # 対数尤度の計算
    LLm <- apply(theta, 1, LL, u=xi, a=a, b=b, c=c)
    Lm <- exp(LLm)
    
    # 事後分布の分母にあたる，定数項
    const <- sum(Lm*w,na.rm = T)
    
    # 事後分布の重み
    Gm <- Lm*w/const
    
    # EAP
    eap <- sum(theta*Gm,na.rm = T)
    
    # 事後標準誤差
    SE <- sqrt(sum((theta-eap)^2*Gm))
    
    # 重みに対応する分点の値をかけて，和を取る＝EAP推定値
    return(c(eap, SE ,const))
  }
  
  if((est == "PVs") || (est == "EAP")){
    # 全受検者のデータをapply関数で与え，EAP推定を実行
    library(pracma)
    eap_apply <- apply(x.all,1,FthetaEAP,a=a,b=b,c=c,theta=xnodes,w=weight)
    message("周辺確率およびEAPの計算が終了しました。")
  }
  
  
  #----------------------------------------------------------#
  #Eatimate MAP
  #----------------------------------------------------------#
  
  FthetaMAP <- function(xi,a,b,c,mu,sigma){
    #----------------------------------------------------------------#
    #MAP推定をおこなう関数。
    #xall     : 項目反応データ
    #a&b      : パラメタファイル。aが識別力，bが困難度
    #mu&sigma : 事前分布の平均と標準偏差
    #----------------------------------------------------------------#
    gid <- xi[1]
    xi <- xi[-1]
    #P(theta) in two-parameter logisticmodel
    ptheta <- function(theta,a,b,c){
      D <- 1.702
      c+(1-c)/(1+exp(-D*a*(theta-b)))
    }
    # 尤度関数
    LL_b <- function(u,theta,a,b,c,mu,sigma){
      sum(log(ptheta(theta,a,b,c))*u+log(1-ptheta(theta,a,b,c))*(1-u)-0.5*((theta-mu)/sigma)^2,na.rm = T)
    }
    
    # 一階偏微分
    fpdLPD <- function(xi, theta,a,b,c,mu,sigma){
      D <- 1.702
      p <- ptheta(theta,a,b,c)
      D*sum(a*(xi - p)*(p-c)/(p*(1-c)), na.rm = T)-1/sigma^2*(theta-mu)
    }
    
    # 二階偏微分
    spdLPD <- function(theta,a,b,c,sigma){
      # 2PLM irf
      D <- 1.702
      p <- ptheta(theta,a,b,c)
      D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T) - 1/sigma
    }
    
    mm <- groupitem[gid]
    
    #初期値を設定。ログオッズ比を用いた。
    if(sum(xi, na.rm = TRUE) == 0){
      t0 <- log(0.5)
    }else if(sum((xi==1)*1,na.rm=T) == length(na.omit(xi))){
      t0 <- log(mm-0.5) 
    }else{
      t0 <- log(sum(xi, na.rm = TRUE)/(mm-sum(xi, na.rm = TRUE)))
    }
    
    d <- 0
    conv1 <- 0
    if(method=="NR"){
      while(conv1==0){
        t1 <- t0 - fpdLPD(xi,t0,a,b,c,mu,sigma)/spdLPD(t0,a,b,c,sigma)
        if(abs(t1-t0)<0.001 || abs(t0)>10 ||is.nan(t1)) {
          conv1 <- 1
        }else{
          t0 <- t1
          d <- d +1
          if(d > 100) { # フィッシャースコアリングの繰り返しが100回を超えたら，二分法に切り替える。
            conv2 <- 0
            p <- maxtheta ; q <- mintheta
            while(conv2 == 0){
              pf <- fpdLPD(xi,p,a,b,c,mu,sigma)
              qf <- fpdLPD(xi,q,a,b,c,mu,sigma)
              t1 <- (p+q)/2
              mf <- fpdLPD(xi,t1,a,b,c,mu,sigma)
              if (abs(pf-qf)<0.001 || c == 1000){
                cat("change bisection method. \n")
                conv2 <- 1
                conv1 <- 1
              }else{
                if(pf*mf>0) {
                  p <- t1
                  d <- d+1
                }else if(qf*mf>0){
                  q <- t1
                  d <- d+1
                } else {
                  conv2 <- 1
                  conv1 <- 1
                }
              }
            }
          }
        }
      }
    } else if(method=="Brent") {
      t1 <- optimise(LL_b,c(mintheta,maxtheta), maximum = T,u=xi,a=a,b=b,c=c,mu=mu,sigma=sigma)$maximum
    } else if(method=="SANN"){
      t1 <- optim(par=t0,fn=LL_b,method="SANN",control=list(fnscale=-1),u=xi,a=a,b=b,c=c,mu=mu,sigma=sigma)$par
    } else {
      t1 <-optim(par=t0,fn=LL_b,gr= function(u, theta,a,b,c,mu,sigma){
        D <- 1.702
        p <- ptheta(theta,a,b,c)
        D*sum(a*(u - p)*(p-c)/(p*(1-c)), na.rm = T)-1/sigma^2*(theta-mu)
      },method=method,control=list(fnscale=-1),u=xi,a=a,b=b,c=c,mu=mu,sigma=sigma)$par
    }
    t1
  }
  
  if((est == "PVs") || (est == "MAP")){
    # 全受検者のデータをapply関数で与え，MAP推定を実行
    map_apply <- apply(cbind(group, x.all),1,FthetaMAP,a=a,b=b,c=c,mu=mu,sigma=sigma)
    message("MAP推定値の計算が終了しました。")
  }
  
  #----------------------------------------------------------#
  #Eatimate MLE
  #----------------------------------------------------------#
  
  FthetaMLE <- function(xi,a,b,c){
    #----------------------------------------------------------------#
    #MLE推定をおこなう関数。
    #xi     : 項目反応データ&母集団変数
    #a&b      : パラメタファイル。aが識別力，bが困難度
    #----------------------------------------------------------------#
    gid <- xi[1]
    xi <- xi[-1]
    #P(theta) in two-parameter logisticmodel
    ptheta <- function(theta,a,b,c){
      D <- 1.702
      c+(1-c)/(1+exp(-D*a*(theta-b)))
    }
    # 対数尤度
    LL <- function(u,theta,a,b,c){
      sum(log(ptheta(theta,a,b,c))*u+log(1-ptheta(theta,a,b,c))*(1-u),na.rm = T)
    }
    # 一階偏微分
    fpd <- function(xi, theta,a,b,c){
      D <- 1.702
      p <- ptheta(theta,a,b,c)
      D*sum(a*(xi - p)*(p-c)/(p*(1-c)), na.rm = T)
    }
    
    # 二階偏微分
    spd <- function(xi,theta,a,b,c){
      D <- 1.702
      p <- ptheta(theta,a,b,c)
      D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T)
    }
    
    d <- 0
    conv1 <- 0
    
    mm <- groupitem[gid]
    temp <- sum(xi==1) == sum(xi,na.rm=T)
    #初期値を設定。ログオッズ比を用いた。
    if(sum(xi, na.rm = T) == 0){
      t1 <- NA
      conv1 <- 1
    }else if(sum((xi==1)*1,na.rm=T) == length(na.omit(xi))){
      t1 <- NA
      conv1 <- 1
    }else{
      t0 <- log(sum(xi, na.rm=T)/(mm-sum(xi, na.rm=T)))
    }
    
    
    if(method=="NR"){
      while(conv1==0){
        t1 <- t0 - fpd(xi,t0,a,b,c)/spd(xi,t0,a,b,c)
        if(abs(t1-t0)<0.001 && !is.nan(t1)) {
          conv1 <- 1
        } else {
          t0 <- t1
          d <- d +1
          if(d > 100 || abs(t0)>10 ||is.nan(t1)){ # ニュートン法の繰り返しが100回を超えたら，二分法に切り替える。
            cat("change bisection method.\n")
            conv2 <- 0
            p <- maxtheta ; q <- mintheta
            while(conv2 == 0){
              pf <- fpd(xi,p,a,b,c)
              qf <- fpd(xi,q,a,b,c)
              t1 <- (p+q)/2
              mf <- fpd(xi,t1,a,b,c)
              if (abs(pf-qf)<0.001 || c == 1000){
                conv2 <- 1
                conv1 <- 1
              }else if(pf*mf>0) {
                p <- t1
                d <- d+1
              }else if(qf*mf>0){
                q <- t1
                d <- d+1
              } else {
                conv2 <- 1
                conv1 <- 1
              }
            }
          }
        }
      }
    } else if(method=="Brent" && conv1 != 1) {
      t1 <- optimise(LL,c(mintheta,maxtheta), maximum = T,u=xi,a=a,b=b,c=c)$maximum
    } else if(method=="SANN" && conv1 != 1){
      t1 <- optim(par=t0,fn=LL,method="SANN",u=xi,control=list(fnscale=-1),a=a,b=b,c=c)$par
    } else if(conv1 != 1){
      t1 <-optim(par=t0,fn=LL,gr=function(u, theta,a,b,c){
        D <- 1.702
        p <- ptheta(theta,a,b,c)
        D*sum(a*(u - p)*(p-c)/(p*(1-c)), na.rm = T)
      },method=method,u=xi,control=list(fnscale=-1),a=a,b=b,c=c)$par
    }
    t1
  }
  
  if(est == "MLE"){
    # 全受検者のデータをapply関数で与え，MAP推定を実行
    mle_apply <- apply(cbind(group, x.all),1,FthetaMLE,a=a,b=b,c=c)
    message("MLE推定値の計算が終了しました。")
  }
  
  
  #----------------------------------------------------------#
  #Eatimate WLE
  #----------------------------------------------------------#
  
  FthetaWLE <- function(xi,a,b,c){
    #----------------------------------------------------------------#
    #重み付け最尤推定をおこなう関数。
    #xall     : 項目反応データ
    #a&b      : パラメタファイル。aが識別力，bが困難度
    #----------------------------------------------------------------#
    
    #P(theta) in two-parameter logisticmodel
    ptheta <- function(theta,a,b,c){
      D <- 1.702
      c+(1-c)/(1+exp(-D*a*(theta-b)))
    }
    
    # テスト情報量（尤度関数の二階偏微分の負の期待値）
    pitheta <- function(theta,a,b,c){
      # 2PLM irf
      ptheta <- function(theta,a,b,c){
        D <- 1.702
        c+(1-c)/(1+exp(-D*a*(theta-b)))
      }
      #
      D <- 1.702
      p <- ptheta(theta,a,b,c)
      I <- D^2*sum(a^2*(1-p)*(p-c)^2/((1-c)^2*p), na.rm = T)
      1/sqrt(I)
    }
    
    # likelihood function with weight
    WL <- function(xi, theta, a, b,c){
      D <- 1.702
      p <- ptheta(theta,a,b,c)
      pi <- pitheta(theta,a,b,c)
      L1 <- sum(xi*log(p)+(xi-1)*log(1-p),na.rm = T)
      L2 <- sum(xi*log(p),na.rm=T)
      B <- sqrt(pi)
      L1 + L2 + B
    }
    
    
    #初期値を設定。ログオッズ比を用いた。
    if(sum(xi, na.rm = TRUE) == 0){
      t0 <- log(0.5)
    }else if (sum(xi, na.rm = TRUE) == m){
      t0 <- log(m-0.5) 
    }else{
      t0 <- log(sum(xi, na.rm = TRUE)/(m-sum(xi, na.rm = TRUE)))
    }
    
    opt <- optimise(WL,interval = c(mintheta,maxtheta),maximum = T,xi=xi,a=a,b=b,c=c)
    t1 <- opt$maximum
    t1
  }
  if(est == "WLE"){
    # 全受検者のデータをapply関数で与え，MAP推定を実行
    wle_apply <- apply(x.all,1,FthetaWLE,a=a,b=b,c=c)
    message("WLEの計算が終了しました。")
  }
  
  
  #----------------------------------------------------------#
  # person fit index: z_3 statistics
  #----------------------------------------------------------#
  pfit <- function(dat, a, b, c){
    #----------------------------------------------------------#
    # 個人適合度を推定するための関数
    # dat : 一列目に最尤推定値，二列目以降に項目反応データを格納した行列，もしくはデータフレーム
    # a$b : 項目パラメタ
    #----------------------------------------------------------#
    #decompose dat
    theta <- dat[1]
    xi <- dat[-1]
    
    # IRF
    ptheta <- function(theta,a,b,c){
      D <- 1.702
      c+(1-c)/(1+exp(-D*a*(theta-b)))
    }
    
    # log likelihood function
    LL <- function(u,theta,a,b,c){
      sum(
        u*log(ptheta(theta,a,b,c))+(1-u)*log(1-ptheta(theta,a,b,c))
        ,na.rm = T
      )
    }
    
    # E_3(theta) : expectation
    E3 <- function(theta, a, b,c){
      sum(
        ptheta(theta, a, b,c) * log(ptheta(theta, a, b,c)) +
          (1 - ptheta(theta, a, b,c)) * log(1- ptheta(theta, a, b,c))
        ,na.rm = T
      )
    }
    
    # sigma3
    S3 <- function(theta, a, b,c){
      sqrt(
        sum(
          ptheta(theta, a, b,c) * (1 - ptheta(theta, a, b,c)) * 
            (log(ptheta(theta, a, b,c)/(1 - ptheta(theta, a, b,c))))^2
          ,na.rm = T
        )
      )
    }
    
    (LL(xi, theta, a, b,c) - E3(theta, a, b,c)) / S3(theta, a, b,c)
  }
  
  lolF <- function(dat,a,b,c){
    theta <- dat[1]
    u <- dat[-1]
    p <- c+(1-c)/(1+exp(-1.702*a*(theta-b)))
    LL <- sum(u*log(p)+(1-u)*log(1-p),na.rm = T)
  }
  
  
  if(est == "EAP"){
    
    # EAP
    eap_m <- eap_apply[1,] %>% round(digits = 5)
    # estimate standard error for EAP estimator
    SE    <- eap_apply[2,] %>% round(digits = 5)
    z3 <- apply(cbind(eap_apply[1,], x.all), 1, pfit, a=a,b=b,c=c) %>% round(digits = 5)
    lol <-  apply(cbind(eap_apply[1,], x.all), 1, lolF, a=a,b=b,c=c) %>% round(digits = 5)
    
    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,EAP=eap_m,SE=SE,z3=z3,lol=lol)
    list("res" = result, "EAPmean&sd" = c(mean(eap_m),sd(eap_m)) %>% round(digits = 5)) 
    
  }else if(est == "MAP"){
    
    #estimate standard error for MAP estimator
    pitheta <- function(dat,a,b,c,sigma){
      # 2PLM irf
      theta <- dat[1]
      xi <- dat[-1]
      ptheta <- function(theta,a,b,c){
        D <- 1.702
        c+(1-c)/(1+exp(-D*a*(theta-b)))
      }
      #
      D <- 1.702
      p <- ptheta(theta,a,b,c)
      I <- -D^2*sum(a^2*(p-c)*(xi*c-p^2)*(1-p)/(p^2*(1-c)^2),na.rm = T) + 1/sigma
      1/sqrt(I)
    }
    SE <- apply(cbind(map_apply, x.all), 1, pitheta, a=a,b=b,c=c,sigma=sigma) %>% round(digits = 5)
    z3 <- apply(cbind(map_apply, x.all), 1, pfit, a=a,b=b,c=c) %>% round(digits = 5)
    lol <- apply(cbind(map_apply, x.all), 1, lolF, a=a,b=b,c=c) %>% round(digits = 5)
    
    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,MAP=map_apply %>% round(digits = 5), SE=SE, z3=z3, lol=lol)
    list("res" = result, "MAPmean&sd" = c(mean(map_apply), sd(map_apply)) %>% round(digits = 5))
    
  }else if(est == "WLE"){
    
    z3 <- apply(cbind(wle_apply, x.all), 1, pfit, a=a,b=b,c=c) %>% round(digits = 5)
    lol <- apply(cbind(wle_apply, x.all), 1, lolF, a=a,b=b,c=c) %>% round(digits = 5)
    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,WLE=wle_apply,z3=z3,lol=lol)
    list("res" = result, "WLEmean&sd" = c(mean(wle_apply), sd(wle_apply)))
    
  }else if(est == "MLE"){
    
    # estimate standard error for MLE estimator
    pitheta <- function(dat,a,b,c){
      theta <- dat[1]
      xi <- dat[-1]
      # 2PLM irf
      ptheta <- function(theta,a,b,c){
        D <- 1.702
        c+(1-c)/(1+exp(-D*a*(theta-b)))
      }
      #
      D <- 1.702
      p <- ptheta(theta,a,b,c)
      I <- D^2*sum(a^2*(1-p)*(p-c)^2/((1-c)^2*p), na.rm = T)
      1/sqrt(I)
    }
    SE <- apply(matrix(mle_apply, ncol = 1), 1, pitheta, a=a, b=b,c=c) %>% round(digits = 5)
    
    # person fit index
    z3 <- apply(cbind(mle_apply, x.all), 1, pfit, a=a,b=b,c=c) %>% round(digits = 5)
    lol <- apply(cbind(mle_apply, x.all), 1, lolF, a=a,b=b,c=c) %>% round(digits = 5)
    
    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,MLE=mle_apply %>% round(digits = 5), SE=SE, z3=z3,lol=lol)
    list("res" = result, "MLEmean&sd" = c(mean(mle_apply, na.rm = T), sd(mle_apply, na.rm = T)) %>% round(digits = 5))
    
    
  }else if( est == "PVs"){
    
    #----------------------------------------------------------#
    # Max. Density of Posterior Distributions
    #----------------------------------------------------------#
    
    Fmaxpdc <- function(xi,theta,a,b,c){
      #----------------------------------------------------------------#
      # 与えられたθに対応する確率密度を計算するための関数。
      # xi    : 受検者iの項目反応データ
      # theta : 潜在特性値
      # a&b   : 項目パラメタ
      #----------------------------------------------------------------#
      theta <- as.matrix(theta)
      xi <- as.numeric(xi)
      ptheta <- function(theta, a, b,c){
        D <- 1.702
        c+(1-c)/(1+exp(-D*a*(theta-b)))
      }
      LL <- function(u,theta,a,b,c){
        sum(u*log(ptheta(theta,a,b,c))+(1-u)*log(1-ptheta(theta,a,b,c)),na.rm = T)
      }
      
      # 項目反応データをapplyで与えて，対数尤度を計算する。
      LLm <- apply(theta,1,LL,u=xi,a=a,b=b,c)
      exp(LLm)
    }
    
    #----------------------------------------------------------#
    #the function for "fg"
    #----------------------------------------------------------#
    
    Ffg <- function(xi,theta,a,b,c){
      #----------------------------------------------------------------#
      # 事後分布の分子にあたる部分を計算するための関数。ただし，事前分布は平均0，標準偏差1
      # xi    : 受検者iの項目反応データ
      # theta : 事前分布
      # a&b   : 項目パラメタ
      #----------------------------------------------------------------#
      exp(sum(xi*log( c+(1-c)/(1+exp(-1.7*a*(theta-b)))) + (1-xi)*log(1- (c+(1-c)/(1+exp(-1.7*a*(theta-b))))), na.rm=T))*dnorm(theta,mean=mu,sd=sigma)
    }
    
    
    #--------------------------------------------------#
    # Rejection Samling starts here.
    #--------------------------------------------------#
    
    # 使用するベクトルと行列を予め作成しておく。
    pv <- matrix(0,n,nofrands) 
    times <- 0
    
    message("フォンノイマンの棄却法にもとづいて推算値のランダムサンプリングを開始します。")
    
    for(k in 1:n){
      #--------------------------------------------------#
      # 推算値算出の流れ
      # 1)乱数発生時の定義域と値域を設定する。設定にはそれぞれEAP推定値とMAP推定値を参考にする。
      # 2)θ軸の乱数(y)とP(θ)軸の乱数(z)を，それぞれ一様分布より得る。
      # 3)項目反応データ，項目パラメタ，yの値を用いて，事後確率密度を計算する(fgvalue)。
      # 4)fgvalueの値よりもzの値が小さければ採択，zの値の方が大きければ棄却する。
      # 5)nofrandsで設定した個数分，推算値を採択できればサブルーチンを終了する。
      # 6)全受検者で1から5を繰り返す。
      #--------------------------------------------------#
      xi <- x.all[k,]
      times <- times + 1
      
      #すでに計算してあるEAP推定値と，事後分布の分子を取り出して行列として展開。
      eap   <- eap_apply[1,k]
      const <- eap_apply[3,k]
      
      #乱数発生時の，P(θ)軸の最大値を設定。
      yheight <- Fmaxpdc(xi,map_apply[k],a,b)*1.001
      
      # 乱数発生時の，θ軸の最大値を設定。
      zmin <- -maxtheta + eap
      zmax <- mintheta + eap
      
      nofpv <- 0
      times_sub <- 0
      while( nofpv <= nofrands ){ 
        
        y <- runif( 1, 0, yheight)
        z <- runif( 1, zmin, zmax)
        times_sub <- times_sub +1
        fg <- apply(xi,1,Ffg,theta=z,a=a,b=b,c=c)
        fgvalue <- fg/const
        
        if( y <= fgvalue){ 
          nofpv <- nofpv + 1
          if( nofpv > nofrands) break
          pv[k,nofpv] <- z
        } 
      } 
      
      if(times==1000){
        message(k,"人目の受検者の推算値の発生が終了しました。")
        times <- 0
      }
      
    } 
    message("推算値の計算が終了しました。")
    
    #推算値の組ごとにもとめた統計量の平均を計算。なお，母集団が複数存在する場合には母集団ごとに求める。
    pvG <- data.frame(group,pv)
    group_mean <- matrix(0,G,nofrands)
    group_var <- matrix(0,G,nofrands)
    group_sd <- matrix(0,G,nofrands)
    
    for (i in 1:G){
      # count the number of subjects for each groups
      ng[i] <- nrow(pvG)
      
      # estiamte group statistics
      group_pv <- pvG %>% dplyr::filter(group==i) %>% dplyr::select(-group)
      group_mean[i,] <- apply(group_pv,2,mean)
      group_var[i,] <- apply(group_pv,2,var)
      group_sd[i,] <- apply(group_pv,2,sd)
    }
    
    M_M <- apply(group_mean,1,mean)
    M_V <- apply(group_var,1,mean)
    M_SD <- apply(group_sd,1,mean) 
    
    PS <- data.frame(group=c(1:G),N=ng,mean=M_M, variance=M_V, sd=M_SD)
    
    message("推算値の集団統計量の計算が終了しました。")
    
    #Little & Rubin. (2002). にもとづく補完間分散
    #V_IMP=(1+1/K)[1/(K-1)sigma_i[(M_PVs-MM_PVs)^2]]+1/K sigma_i[V(M_PV_i)]
    
    MI_SE <- function(M){#var() is the function for unbias variance
      K <- length(M)
      (K+1)/K*var(M) + (K-1)/K*var(M)
    }
    
    SE_M <- apply(group_mean,1,MI_SE)
    SE_V <- apply(group_var,1,MI_SE)
    SE_SD <- apply(group_sd,1,MI_SE)
    
    SE <- data.frame(group=c(1:G),mean=SE_M, variance=SE_V, sd=SE_SD)
    
    message("集団統計量の標準誤差の推定が終了しました。")
    
    #推算値の平均(Right & Wrong)
    pvmeans <- pvmeans_w <- apply(pv,1,mean)
    pvmeans_r <- apply(pv,2,mean)
    
    #result
    
    result <- data.frame(ID=ID,GROUP=group,SCORE=xscore,EAP=eap_apply[1,],MAP=map_apply,PVmeans_W=pvmeans_w,PV=pv,AREAS=const)
    plausible_values <-data.frame(ID=ID,Group=group,PV=pv)
    if(output==T){
      message("結果のCSVファイルを書き出しています。")
      write.csv(result,paste0(file,"_result.csv"),quote=F,row.names = F)
      write.csv(plausible_values,paste0(file,"_PVs.csv"),quote=F,row.names = F)
      write.csv(PS,paste0(file,"_PVS population statistics.csv"),quote=F,row.names = F)
      write.csv(SE,paste0(file,"_PVS standard error.csv"),quote=F,row.names = F)
    }
    pv <- data.frame(ID=ID,group=group,SCORE=xscore,EAP=eap_apply[1,],MAP=map_apply,PV=pv)
    list("PVs" = pv, "EAPmean&sd" = c(mean(eap),sd(eap)), 
         "MAPmean&sd" = c(mean(map_apply), sd(map_apply)), "PVmean&sd" = c(M_M, M_SD), "PS" = PS, "SE" = SE) 
  }
}



#-------------------------------------------------------------------------------------------------------#
#シミュレーションデータ生成関数
#-------------------------------------------------------------------------------------------------------#


sim_gen <- function(theta, a, b, item = 'A', root=1/5){
  
  #　あらかじめ能力値とパラメタは乱数ベクトルを用意しておく。
  ptheta <- function(theta, a, b){
    # IRT 2 PLM response probability
    D <- 1.702
    1/(1+exp(-D*a*(theta-b)))
  }
  # ptheta(theta[1], a, b)
  prob <- t(apply(matrix(theta, ncol = 1), 1, ptheta, a =a, b = b))
  resfunc <- function(prob,root=root){
    # 反応確率と一様乱数から01データを発生させる関数。
    # 受検者一人分の正答確率を与える。（apply関数などで）
    # rootの値をいじることで，NAの発生率を変更できる。rootの値が小さいほど，発生率も小さい。
    prob <- matrix(prob, ncol = 1)
    subfunc <- function(prob){
      if(prob < runif(1)) res <- 0
      else res <- 1
      sample(x=c(res,NA),size=1,prob=c(prob^root,1-prob^root))
    }
    res <- apply(prob, 1, subfunc)
    return(res)
  }
  
  # generate item respinse pattern
  test <- t(apply(prob, 1, resfunc))
  # change item name
  itemn <- formatC(c(1:length(a)), width = 3, flag = 0)
  itemn <- apply(matrix(item, ncol = 1), 1, paste0, itemn)
  colnames(test) <- itemn
  return(test)
}





#-------------------------------------------------------------------------------------------------#
#等化係数推定関数
#-------------------------------------------------------------------------------------------------#

CEquating <- function(T_para, F_para, method="HB", D =1.702, graph = F, Fname ="default", output = F, Change = 0){
  #--------------------------------------#
  #parameter file format is adapted to EasyEstimation Para.csv
  #When you read parameter file, you should specify argument"skip = 1" and "header = F".
  #the implement of graph function has been incompleted now.Don't use.
  #--------------------------------------#
  
  Tpara <- T_para[T_para$V3 != 0,]   #remove item that a-parameter is 0
  Fpara <- F_para[F_para$V3 != 0,]
  CII <- Fpara[Fpara$V1 %in% Tpara$V1,1]  #Comon Item Index
  at <- bt <- af <- bf <- numeric(length(CII))
  
  for (i in 1:length(CII)){
    at[i] <- Tpara[Tpara$V1 == CII[i], "V3"]
    bt[i] <- Tpara[Tpara$V1 == CII[i], "V4"]
    af[i] <- Fpara[Fpara$V1 == CII[i], "V3"]
    bf[i] <- Fpara[Fpara$V1 == CII[i], "V4"]
  }
  
  #mean & sigma method
  tm <- mean(bt)  #mean of diff para
  ts <- sqrt(mean((bt-tm)^2))    #standard deviation(It's not used unbiased estimator)
  fm <- mean(bf)
  fs <- sqrt(mean((bf-fm)^2))
  msA <- ts/fs
  msK <- tm - fm * msA
  
  p <- seq(-4,4,8/30)   #nodes of division quadrature
  w <- numeric(31)      #weights of division quadrature
  for (m in 1:31){    w[m] <- dnorm(p[m],0,1)*(8/30)  }
  
  if(method == "HB"){
    #criterion function
    cf <- function(x,af,bf,at,bt,D=D){
      A <- x[1]
      K <- x[2]
      Q1 <- 0
      Q2 <- 0
      ptheta <- function(theta, a, b, D){
        1/(1+exp(-D*a*(theta-b)))
      }
      for(m in 1:31){
        HC1 <- (ptheta(p[m],at,bt,D)-ptheta(p[m],af,bf,D))^2
        HC2 <- (ptheta(p[m],af,bf,D)-ptheta(p[m],at*A,(bt-K)/A,D))^2
        
        Q1 <- Q1 + sum(HC1*w[m])
        Q2 <- Q2 + sum(HC2*w[m])
      }
      Q <- Q1 +Q2
      return(Q)
    }
    res <- nlm(f=cf, p=c(msA,msK),af=af,bf=bf,at=at,bt=bt,D=D)#initial value is equating coefficient estimated by mean & sigma
    
    A <- res$estimate[1]
    K <- res$estimate[2]
    
      if(graph == T){
        n <- ceiling(sqrt(length(CII)))
        par(mfrow=c(n,n))
        for(j in length(CII)){
          ptheta <- function(theta, a, b, D){
            1/(1+exp(-D*a*(theta-b)))
          }
          curve(ptheta(x, at[j], bt[j],D),-4,4,ylim=c(0:1),lty=1,lwd=2,main=CII[j])
          par(new=T)
          curve(ptheta(x, af[j]/A, (A*bf[j]+K),D),-4,4,ylim=c(0:1),lty=2,lwd=2,main="")
          par(new=T)
          curve(ptheta(x, af[j], bf[j],D),-4,4,ylim=c(0:1),lty=3,lwd=2,main="")
          legend(-4,1.0,c("T","F*","F"),lty=1:3,lwd=2)
        }
      }
  } else if (method == "SL"){
    cf2 <- function(x,af,bf,at,bt,D){
      A <- x[1]
      K <- x[2]
      Q1 <- 0
      Q2 <- 0
      ptheta <- function(theta, a, b, D){
        1/(1+exp(-D*a*(theta-b)))
      }
      for(m in 1:31){
        SLC1 <- (sum(ptheta(p[m],at,bt,D))-sum(ptheta(p[m],af,bf,D)))^2
        SLC2 <- (sum(ptheta(p[m],af,bf,D))-sum(ptheta(p[m],at*A,(bt-K)/A,D)))^2
        
        Q1 <- Q1 + SLC1*w[m]
        Q2 <- Q2 + SLC2*w[m]
        
      }
      Q <- Q1 +Q2
      return(Q)
    }
    res <- nlm(f=cf2, p=c(msA,msK),af=af,bf=bf,at=at,bt=bt,D=D)
    
    A <- res$estimate[1]
    K <- res$estimate[2]
    
      if(graph == T){#incompleted option
        tcc <- function(theta, a, b, D){
          res <- numeric(length(theta))
          for(i in 1:length(theta)){
            res[i] <- sum(ptheta(theta[i], a, b, D))
          }
          return(res)
          
      }
    }
  }
  
  #equating Fpara on the scale of Tpara
  nCII <- numeric(0)
  if(nrow(Fpara) != length(CII)){
    nCII <- Fpara[!(Fpara$V1 %in% Tpara$V1),1]  #extract non comon item 
    #common items
    for (n in nCII){
      Fpara[Fpara$V1==n,"V3"] <- round(Fpara[Fpara$V1==n,"V3"]/A, digits = 5)
      Fpara[Fpara$V1==n,"V4"] <- round(Fpara[Fpara$V1==n,"V4"]*A+K, digits = 5)
    }
  }
  #non comon items
  for (i in CII){
    if(Change == 0){
      Fpara[Fpara$V1==i,"V3"] <- Tpara[Tpara$V1==i,"V3"]
      Fpara[Fpara$V1==i,"V4"] <- Tpara[Tpara$V1==i,"V4"]
    } else if(Change == 1){
      Fpara[Fpara$V1 == i,"V3"] <- round(sqrt(Fpara[Fpara$V1 == i, "V3"]/A * Tpara[Tpara$V1 == i, "V3"]), digits = 5) # geometric mean
      Fpara[Fpara$V1 == i,"V4"] <- round((Fpara[Fpara$V1 == i, "V4"]*A+K + Tpara[Tpara$V1 == i, "V4"])/2, digits = 5) # arithmetic mean
    } else if(Change == 2){
      Fpara[Fpara$V1 == i,"V3"] <- round(Fpara[Fpara$V1 == i, "V3"]/A, digits = 5) # geometric mean
      Fpara[Fpara$V1 == i,"V4"] <- round(Fpara[Fpara$V1 == i, "V4"]*A+K, digits = 5) # arithmetic mean
    }
  }
  
  #replace T_para to Tpara exept the a=0 item
  #for (j in Fpara$V1){
  #  F_para[F_para$V1 == j,"V3"] <- Fpara[Fpara$V1 == j, "V3"]
  #  F_para[F_para$V1 == j,"V4"] <- Fpara[Fpara$V1 == j, "V4"]
  #}

  result <- list(D=paste0("D = ",D), para=Fpara, METHOD=method, 
                 EquatingCoefficient_A_K=c(A,K), MeanSigma_A_K=c(msA,msK),ComonItem=CII,NonComonItem=nCII)
  
  if(output == T){
    write.table(result[1],file=paste0(Fname,"ParaEquated.csv"),sep = "",quote = F, col.names = F, row.names = F)
    write.table(result[2],file=paste0(Fname,"ParaEquated.csv"),sep = ",",quote = F, append = T, col.names = F, row.names = F)
  }

  
  return(result)
  
}


#-------------------------------------------------------------------------------------------------#
# Library function
#-------------------------------------------------------------------------------------------------#


Library <- function(x){ # 引数は""で囲むこと
  tryCatch(
    library(x, character.only = T),
    error = function(x){
      install.packages(x)
      library(x, character.only = T)
    }, 
    finally = message(x,"package installing is success")
  )
}



#-------------------------------------------------------------------------------------------------#
# timer function
#-------------------------------------------------------------------------------------------------#

timer <- function(StopTime, MESSAGE = "STOP!!", disp = T, units = "secs", 
                  cex = 5, # だいたい10位がベスト 
                  font = 1, 
                  col = 1, 
                  family="",  # Macの場合は，Helveticaを選択すればOK，ただし日本語は不可。
                  os = "win"
){
  time <- 0
  if(disp == T){
    repeat {
      if(time == 0) START <- Sys.time()
      t0 <- Sys.time()
      DIFF <- 0
      while(DIFF == 0){
        t1 <- Sys.time()
        DIFF <- t1 - t0
        if(DIFF < 0.01) DIFF <- 0
      }
      time <- t1 - START
      if(time >= StopTime) break
      time <- round(as.numeric(time, units = units), digits = 3)
      cat(time, units,".\r")
    }
    
    cat("\n",MESSAGE)
    
    if(os =="win")    win.graph()
    if(os == "mac") quartz()
    
    plot(NA, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    par(cex =cex, col = col, font = font, family = family)
    text(0.5,0.5,MESSAGE)
  }
}




#-------------------------------------------------------------------------------------------------#
# 被検者対数尤度計算関数
#-------------------------------------------------------------------------------------------------#

Flol <- function(theta,x,a,b,c){
  #cat("識別力が0の項目を削除します。\n")
  x <- x[,a!= 0]
  #cat(sum((a == 0)*1),"個の項目が削除されました。\n")
  c <- c[a!=0]
  b <- b[a!=0]
  a <- a[a!=0]
  dat <- cbind(theta,x)
  lolF <- function(dat,a,b){
    theta <- dat[1]
    u <- dat[-1]
    p <- c+(1-c)/(1+exp(-1.702*a*(theta-b)))
    LL <- sum(u*log(p)+(1-u)*log(1-p),na.rm = T)
  }
  apply(dat, 1, lolF, a=a, b=b)
}



#-------------------------------------------------------------------------------------------------#
# 復元得点分布算関数
#-------------------------------------------------------------------------------------------------#
obscore_dist <- function(theta,a,b,name="test",color="cyan", output=1){ 
  #--------------------------------------------------------------#
  # This function generates distribution of IRT observed scores
  # on 2-parameter logistic model.
  #   
  #　　m　　　	: 　 number of items
  #　　n　　　　:  　number of subjects　
  #　　rn          :　　 number of normal random numbers
  #　　theta	:  　ability
  #　　a　　　　: 　 item discriminating powers
  #　　b　　　　:  　item difficulties
  #　　name    :　　main title of the histogram(base::hist())
  #　　color   :　　color of the histogram
  #    output  :    if 1 return the vector of observed score. if 2 return the histgram.
  #
  #--------------------------------------------------------------#
  
  a <- as.matrix(a)　　　　　　　　　　　　　　　　　#項目パラメタを行列形式にする
  b <- as.matrix(b)　　　
  theta <- as.matrix(theta)
  
  
  #--------------------------------------------------------------#
  # Number of Items and Subjects
  #--------------------------------------------------------------#
  m <- nrow(a)	    #n of items 
  m1 <- m+1
  
  n <- nrow(theta)　	#n of subjects (thetaに推定した能力値を代入する場合)
  
  
  #--------------------------------------------------------------#
  # Distribution of Ability 
  #--------------------------------------------------------------#	
  i <- sort.list(theta)　　　　　　#thetaを昇順に並べる
  theta <- theta[i]
  i <- matrix(1:n,1,n)
  
  
  #--------------------------------------------------------------#
  # Conditional Probabilities of Test Scores
  #--------------------------------------------------------------#
  
  probability <- function(trait,a,b){	
    #--------------------------------------------------------------#
    # This function generates conditional probability
    # for a fixed ability on 2-parameter logistic model.
    #
    #       The original Fortran77 program was 
    #               developed by Inoue,S. , December 1990.
    #               extended by Shibayama,T., January 1991.
    #               translated into R by Shibayama,T. September 2008.
    #               functionalized by Itamiya, C., & Shibuya. T., June 2018.
    #
    #	　m	　:  number of items
    #	　trait	　:  ability
    #	　a	　:  item discriminating powers
    #	　b	　:  item difficulties
    #	　prb	:  conditional probability for a fixed ability
    #	　ptheta	　:  P(theta)
    #	　qtheta	　:  Q(theta)=1-P(theta)
    #--------------------------------------------------------------#
    
    m <- nrow(a)
    m1 <- m+1
    ptheta <- matrix(0,m,1)
    qtheta <- matrix(0,m,1)
    prb <- matrix(0,m1,1)
    #
    ptheta<- 1/(1+exp(-1.7*a*(trait-b)))
    qtheta<- 1-ptheta
    #
    prb[1] <- qtheta[1]
    prb[2] <- ptheta[1]
    #
    for(j in 2:m){
      l <- j -1
      j1 <- j+1
      l1 <- l+1
      prb[j1] <- prb[l1]*ptheta[j]
      
      for(i in l:1){
        k <- i -1
        i1 <- i+1
        k1 <- k+1
        prb[i1]<-prb[k1]*ptheta[j]+prb[i1]*qtheta[j]
      }
      
      prb[1] <- prb[1]*qtheta[j]
    }
    #
    probability <- prb
    #
  }
  
  prbtestscore<-matrix(0,n,m1)
  for(i in 1:n){
    prbtestscore[i,]<- t(probability(theta[i],a,b))
  }
  
  
  #--------------------------------------------------------------#
  # Marginal Probabilities of Test Score
  #--------------------------------------------------------------#
  freq <- t(prbtestscore) %*% matrix(1,n,1)
  freq <- cbind(matrix(0:m,m1,1),freq)　　　　　#[得点｜度数分布]
  temp <- round(freq[,2])					#度数の小数点以下を四捨五入する
  score <- rep(freq[,1],temp)				#度数分布表をベクトルに展開する
  
  if(output == 1){
    return(score)
  }else if(output==2){
    mx <- max(score)		
    hist(score,freq=FALSE,ylim=c(0,0.15),breaks=seq(-0.5,(mx+0.5),1),col=color,main=name,cex.main=1.5)   #相対度数分布のヒストグラムを近似的に描画する
  }
}

#----------------------------------------------------------------#
# true score frequency
#----------------------------------------------------------------#
tscore_dist <- function(theta,a,b){
  #--------------------------------------------------------------#
  a <- as.matrix(a)　　　　　　　　　　　　　　　　　#項目パラメタを行列形式にする
  b <- as.matrix(b)　　　
  theta <- as.matrix(theta) 
  
  # Number of Items and Subjects
  m <- length(a) #n of items
  m1 <- m+1
  n <- length(theta) #n of subjects
  
  
  # Distribution of Ability
  i <- sort.list(theta)
  theta <- theta[i]
  i <- matrix(1:n,1,n)
  
  # Comditional Probabilities of True Scores
  #------------------------------------------------------------#
  truescore<-function(trait,a,b){
    #----------------------------------------------------------#
    # This function generates true scores
    # 2-parameter logistic model.
    #
    ##  m   :  number of items
    #   trait   :  ability
    #   a   :  item discriminating powers
    #   b   :  item difficulties
    #   ptheta  :  P(theta)
    #   tscore  :  true scores
    #
    m <- nrow(a)
    m1 <- m+1
    #ptheta<- matrix(0,m,1)
    #
    ptheta<- 1/(1+exp(-1.7*a*(trait-b)))
    #
    truescore <- sum(ptheta)
    #
  }
  #------------------------------------------------------------#
  tscore <- matrix(0,n,1)
  
  for(i in 1:n){
    tscore[i,] <- t(truescore(theta[i],a,b))
  }
  
  return(round(tscore))
}


#--------------------------------------------------------------#
obscore_dist_data <- function(theta,a,b){ 
  #--------------------------------------------------------------#
  # This function generates distribution of IRT observed scores
  # on 2-parameter logistic model.
  #   
  #m	:  number of items
  #n:  number of subjects
  #rn          : number of normal random numbers
  #theta	:  ability
  #a:  item discriminating powers
  #b:  item difficulties
  #name    :main title of the histogram
  #color   :color of the histogram
  #
  #--------------------------------------------------------------#
  
  a <- as.matrix(a)#鬆?逶ｮ繝代Λ繝｡繧ｿ繧定｡悟?怜ｽ｢蠑上↓縺吶ｋ
  b <- as.matrix(b)
  theta <- as.matrix(theta)
  
  
  #--------------------------------------------------------------#
  # Number of Items and Subjects
  #--------------------------------------------------------------#
  m <- nrow(a)	    #n of items 
  m1 <- m+1
  
  n <- nrow(theta)	#n of subjects (theta縺ｫ謗ｨ螳壹＠縺溯?ｽ蜉帛､繧剃ｻ｣蜈･縺吶ｋ蝣ｴ蜷?)
  
  
  #--------------------------------------------------------------#
  # Distribution of Ability 
  #--------------------------------------------------------------#	
  i <- sort.list(theta)#theta繧呈??鬆?縺ｫ荳ｦ縺ｹ繧?
  theta <- theta[i]
  i <- matrix(1:n,1,n)
  
  
  #--------------------------------------------------------------#
  # Conditional Probabilities of Test Scores
  #--------------------------------------------------------------#
  
  probability <- function(trait,a,b){	
    #--------------------------------------------------------------#
    # This function generates conditional probability
    # for a fixed ability on 2-parameter logistic model.
    #
    #       The original Fortran77 program was 
    #               developed by Inoue,S. , December 1990.
    #               extended by Shibayama,T., January 1991.
    #               translated into R by Shibayama,T. September 2008.
    #
    #	m	:  number of items
    #	trait	:  ability
    #	a	:  item discriminating powers
    #	b	:  item difficulties
    #	prb	:  conditional probability for a fixed ability
    #	ptheta	:  P(theta)
    #	qtheta	:  Q(theta)=1-P(theta)
    #--------------------------------------------------------------#
    
    m <- nrow(a)
    m1 <- m+1
    ptheta <- matrix(0,m,1)
    qtheta <- matrix(0,m,1)
    prb <- matrix(0,m1,1)
    #
    ptheta<- 1/(1+exp(-1.7*a*(trait-b)))
    qtheta<- 1-ptheta
    #
    prb[1] <- qtheta[1]
    prb[2] <- ptheta[1]
    #
    for(j in 2:m){
      l <- j -1
      j1 <- j+1
      l1 <- l+1
      prb[j1] <- prb[l1]*ptheta[j]
      
      for(i in l:1){
        k <- i -1
        i1 <- i+1
        k1 <- k+1
        prb[i1]<-prb[k1]*ptheta[j]+prb[i1]*qtheta[j]
      }
      
      prb[1] <- prb[1]*qtheta[j]
    }
    #
    probability <- prb
    #
  }
  
  prbtestscore<-matrix(0,n,m1)
  for(i in 1:n){
    prbtestscore[i,]<- t(probability(theta[i],a,b))
  }
  
  
  #--------------------------------------------------------------#
  # Marginal Probabilities of Test Score
  #--------------------------------------------------------------#
  freq <- t(prbtestscore) %*% matrix(1,n,1)
  freq <- cbind(matrix(0:m,m1,1),freq)#[蠕礼せ?ｽ懷ｺｦ謨ｰ蛻?蟶ゾ
  #sum(freq[,2])
  
  temp <- round(freq[,2])					#度数の小数点以下を四捨五入する
  score <- rep(freq[,1],temp)				#度数分布表をベクトルに展開する
  
  dist <- freq[,2]/n    #逶ｸ蟇ｾ蠎ｦ謨ｰ蛻?蟶?
  #distribution <- cbind(xfreq,xdist)
  
  cumulative <- cumsum(freq[,2])      #邏ｯ遨榊ｺｦ謨ｰ蛻?蟶?
  cum.relative <- cumulative/n#邏ｯ遨咲嶌蟇ｾ蠎ｦ謨ｰ蛻?蟶?
  
  
  distribution <- data.frame(score=freq[,1],frequency=round(freq[,2]),relative.f=dist,cumulative.f=cumulative,cum.relative.f=cum.relative)
  return(distribution)
  #  return(score)
}

#-------------------------------------------------------------------------------------------------#
# 等パーセンタイル等化関数
#-------------------------------------------------------------------------------------------------#
epe <- function(x, y){
  
  #----------------------------------------------------------------------------------
  #
  #		This function generates the Form Y equipercentle equivalent of score x on FormX, eY(x).
  #		eY(x) is called equating function.
  #
  #		x : 			test score of Form X
  #		nx:			number of items on Form X
  #		y : 			test score of Form Y
  #		ny			number of items on Form Y
  #		dist.f:			function for generating the frequency distribution
  #		prfx:			function for generating the percentile rank on score x in test form X
  #		prfy:			function for generating the percentile rank on score y in test form Y
  #		pfL:			the inverse of the percentile rank function.	
  #		pfU:			the inverse of the percentile rank function.
  #		eYx:			the function of the equipercentile equivalent of score x on the scale of Form Y.
  #		eXy:			the function of the equipercentile equivalent of score y on the scale of Form X.
  #		tablex:			
  #		tabley:			
  #		resultx:		
  #		resulty:		
  #
  #
  #----------------------------------------------------------------------------------	
  
  tablex <- x;tabley <- y
  
  prfx <- function(q){
    x <- q
    x1 <- q+1
    rx <- round(x1)
    ddf <- tablex[,3]
    cddf <- tablex[,5]
    prf <- (cddf[rx]+(x-rx+0.5)*(cddf[rx]-cddf[rx-1]))*100
  }
  
  prfy <- function(r){
    x <- r
    x1 <- r+1
    rx <- round(x1)
    ddf <- tabley[,3]
    cddf <- tabley[,5]
    prf <- (cddf[rx]+(x-rx+0.5)*(cddf[rx]-cddf[rx-1]))*100
  }
  
  #n of item on Form X and Y
  
  nx <- max(x$score)
  ny <- max(y$score)
  
  resultx <- matrix(0, nrow(x), 1)
  resulty <- matrix(0, nrow(y), 1)
  resultx[1,1] <- x[1,3] 
  resulty[1,1] <- y[1,3]
  
  z<- 0
  for (i in 1:nx){
    z <- z +1
    resultx[(z+1), 1] <- prfx(z)
  }
  
  z<- 0
  for (i in 1:ny){
    z <- z +1
    resulty[(z+1), 1] <- prfy(z)
  }
  
  #----------------------------------------------------------------------------------
  # equating function
  # test Form Y equipercentile equivalent of score x on Form X
  #----------------------------------------------------------------------------------
  
  tabley2 <- matrix(0, nrow(y), 3)
  tabley2[,1] <- y$score
  tabley2[,2] <- y[,5]      #relative.p
  tabley2[,3] <- resulty    #percentile rank
  tabley2 <- as.data.frame(tabley2)
  
  #	paste(tabley2)
  
  #----------------------------------------------------------------------------------
  #^ "tabley2 " is used for "pf "
  #Dont forget to make "table2" as data frame.
  #
  #define the percentile function that contains yL
  #(This value is smallest value corresponds with above cumulative percent "p"%
  #----------------------------------------------------------------------------------
  
  pfU <- function(p){
    if(p/100 > tabley2$V2[1]){
      a <- tabley2[tabley2$V2 > (p/100),]
      yU <- a[1,1]
      FyU <- a$V2[1]
      b <- tabley2[tabley2$V2 < (p/100),]
      FyU1 <- b$V2[yU]				#yU-1 is not proper value because table b couont "0".
      pfU <- (p/100 - FyU1)/(FyU-FyU1) + yU - 0.5
    }else{
      pfU <- tabley2$V2[1]
    }
    
  }
  
  #----------------------------------------------------------------------------------
  #Next, define the percentile function that contains yL
  #(This value is largest value corresponds with below cumulative percent "p"%)
  #----------------------------------------------------------------------------------
  
  pfL <- function(p){
    if(p/100 > tabley2$V2[1]){
      a <- tabley2[tabley2$V2 > (p/100),]
      yL1 <- a[1,1]
      FyL1 <- a$V2[1]
      b <- tabley2[tabley2$V2 < (p/100),]
      FyL <- b$V2[yL1]
      pfL <-(p/100 - FyL)/(FyL1-FyL) + (yL1 - 1) + 0.5
    }else{
      pfL <- tabley2$V2[1]
    }
    
  }
  
  eYx <- matrix(0, nrow(x), 2)
  for (i in 1:nrow(x)){
    eYx[i,1] <- pfU(resultx[i])
    eYx[i,2] <- pfL(resultx[i])
  }
  
  result <- cbind(x$score, eYx)
  colnames(result) <- c("X raw score", "pfU", "pfL")
  result
}


