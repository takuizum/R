d1 <- c("good","bad","good","good")
d2 <- c("exciting", "exciting")
d3 <- c("good","good","exciting","boring")
d4 <- c("bad","boring","boring","boring")
d5 <- c("bad","exciting","bad")
d6 <- c("bad","bad","boring","boring")

# 各事例ベクトルをリストに内包する。
P <- list("d1"=d1,"d2"=d2,"d3"=d3)
N <- list("d4"=d4,"d5"=d5,"d6"=d6)

D <- list("P"=P, "N"=N)

# 文字列のマッチングにはgrep関数が使える。
grep("good",P)
grep("bad",P)

# Nさんのひとつめの文書についてマッチング
grep("bad",N$d4)
grep("boring",N$d4)


# ただし，この関数の返値は，マッチした要素の位置を示しているに過ぎない。
# 出現数であれば，この返値の長さを数えれば良い。

length(grep("good",P))
length(grep("bad",P))

length(grep("boring",N$d4))

# これでもよい

sum(N$d4 == "boring")

# ただしこれはダメ

N == "bad"




# 対数線形モデル



words <- c("good", "bad", "exciting", "boring")

# 自動的に文書内の重複なしの単語を取り出すこともできる。
words <- unique(unlist(D))


f <- function(d, y, w, learn = 0.1, conv =0.001){
  
  # d : 文書
  # y : クラス数
  # w : 単語
  
  # マッチング関数
  m <- function(d, w){
    
    # 文書内に指定した単語が出現しているかどうかをT or Fで返して，それを01に書き換える。
    
    a <- d == w
    
    if(sum(a) != 0 ){
      a <- 1
    } else {
      a <- 0
    }
    return(a)
  }
  
  
  # すべての文書ベクトルを，ひとつのリストに格納する。
  dd <- unlist(d, recursive = F) # recursive = F　でリストの構造を残したまま，unlistできる。
  
  # 結果代入用の行列を作成
  # そのためにまずは全，文書数をカウントする
  
  v <- length(dd)
  
  mat <- matrix(0, nrow = v * length(d), ncol = length(w) * length(d))
  
  cc <- 0
  
  for(i in 1:y){
    
    n <- length(dd)
    
    # マッチングした値を代入する行数
    if(i == 1){
      rc1 <- 1
      rc2 <- i * length(w)
    } else {
      rc1 <- rc2 + 1
      rc2 <- i * length(w)
    }
    
    for(j in 1:n){
      
      # マッチングした値を代入する列数
      cc <- cc + 1
      
      w <- as.matrix(w)
      
      mat[cc, rc1:rc2] <- apply(w, 1, m, d = dd[[j]])
      
      # matの構造について
      # クラス>文書
      
    }
    
  }
  
  
  
  mat2 <- matrix(0, nrow = v, ncol = length(w) * length(d))
  
  cc <- 0
  
  for(i in 1:y){
    
    n <- length(d[[i]])
    
    # マッチングした値を代入する行数
    if(i == 1){
      rc1 <- 1
      rc2 <- i * length(w)
    } else {
      rc1 <- rc2 + 1
      rc2 <- i * length(w)
    }
    
    for(j in 1:n){
      
      # マッチングした値を代入する列数
      cc <- cc + 1
      
      w <- as.matrix(w)
      
      mat2[cc, rc1:rc2] <- apply(w, 1, m, d = d[[i]][[j]])
      
      # mat2の構造について
      # クラス>文書
      
    }
    
  }
  
  
  
  # 尤度関数の一階偏微分をもちいて，最急勾配法を実行する。
  
  
  # 関数定義
  
  
  dL <- function(dat, w, dat2, C = 1){
    
    # dat  : 各クラス，全部の事例における素性ベクトルを行に持つ，行列。
    # w    : パラメタ
    # dat2 : 特定のクラスの，特定の文書の素性ベクトルを行に持つ，行列
    # C    : 正の定数
    
    # まずは正則化項は無視する。
    
    sec <- function(dat, w){
      
      Zdw <- sum(exp(dat %*% w))
      
      sigmay <- 
        colSums(
          t(
            apply(dat, 1, 
                  function(x){
                    P <- as.vector(exp(x %*% w) * 1/Zdw)
                    as.vector(x) * P
                  }
            )
          )
        )
      
      t(sigmay)
      
    }

    
    
    # 計算結果代入用ベクトル
    
    res <- numeric(0)
    
    
    # Sigma_yの計算
    for(i in 1:nrow(dat2)){
      
      key <- seq(1, nrow(dat), nrow(dat2)) 
      key <- key + i - 1
      
      first <- dat2[i, ]
      
      second <-sec(dat[key,], w)
      
      res <- rbind(res, (first - second))
      
      if(i == nrow(dat2)){
        res <- colSums(res)
      }
      
    }
    
    res <- as.matrix(res)
    
    # 正則化項の計算
    
    reg <- C * w
    
    return(res-reg)
    
    
  }
  
  
  
  
  # 初期値
  w0 <- matrix(0, nrow = ncol(mat), ncol =1) # initial value
  
  # 繰り返し用
  e <- 0
  
  # カウント用
  count <- 0
  
  while(e == 0){
    
    # 学習率は0.1に設定
    w1 <- w0 + learn * dL(dat = mat, w = w0, dat2 = mat2, C = 1)
    
    if(max(abs(w1 - w0)) < conv) e <- 1
    count <- count + 1
    w0 <- w1
    
  }
  
  res <- list("para" = w1, "iteration" = count)
  
  return(res)
  
}



res1 <- f(D, 2, words)




# 学習したパラメタをもちいて文書を分類する。

# 分類したい文書
dn <- c("exciting", "boring")

# 分類するためのクラス
cla <- c("P","N")


classify <- function(para, d, c, w){
  
  # para : 対数線形モデルにより推定したパラメタ
  #    d : 分類したい文書
  #    c : 分類先のクラス（順番に注意）
  #    w : 単語
  
  
  # マッチング関数
  m <- function(d, w){
    
    # 文書内に指定した単語が出現しているかどうかをT or Fで返して，それを01に書き換える。
    
    a <- d == w
    
    if(sum(a) != 0 ){
      a <- 1
    } else {
      a <- 0
    }
    return(a)
  }
  
  
  # 結果代入用の行列を作成
  mat <- matrix(0, nrow = length(c), ncol = nrow(para))
  
  
  for(i in 1:length(c)){
    
    # マッチングした値を代入する行数
    if(i == 1){
      rc1 <- 1
      rc2 <- i * length(w)
    } else {
      rc1 <- rc2 + 1
      rc2 <- i * length(w)
    }
    
    # マッチングした値を代入する列数はiでよい。
    
    mat[i, rc1:rc2] <- apply(as.matrix(w), 1, m, d = d)
    
  }
  
  
  # 内積を計算し，分類する。
  
  res <- mat %*% para
  rownames(res) <- c

  key <- res == max(res)
  
  message("この文書のクラスは", c[key], "です")
  
  return(res)
  
}


classify(res1$para, dn, c("P", "N"), words)
