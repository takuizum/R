# practice for using "mirt" package

install.packages("mirt")

library(mirt)
library(ltm)
library(irtoys)
library(Rcpp)

ngaku16 <- read.csv("C:/Users/sep10_z1vk0al/OneDrive/Documents/12_R/practice data/ngakudata/ngaku16.csv", header = F)
row.names(ngaku16) <- ngaku16[,1]
ngaku16 <- ngaku16[,-1]


# いろいろな関数の確認

expand.table(LSAT7)

# この関数は，テストデータ("LSAT7")のデータから項目反応データを引っ張ってくる関数。
# LSATのテストデータには，項目反応パタンと，そのパタンの受検者の度数が入っている。それを分析可能なデータとして展開するための関数。


# unidimensional IRT
mod1 <- mirt(RESP[,c(-1,-2)], model = 1, itemtype = "2PL", method = "EM", optimizer = "NR", SE.type = 'complete', SE = T)
summary(mod1, simplify=T)

# UIRTモデルもパラメタがほしいときは引数に"IRTpars=TRUE"をいれてやればよい。
coef(mod1)

models <- "F1 = 1-10"
group <- c(rep("G1",1000), rep("G2",1000), rep("G3",1000), rep("G4",1000), rep("G5",1000))
mirt(RESP[,c(-1,-2)], model = models, itemtype  = "2PL")
res_multi <- multipleGroup(RESP[,c(21:30+2)], models, group = group)
# 集団内に完全に欠測を含んだ項目があると，計算できない模様
coef(res_multi,IRTpars=TRUE)

DIF(res_multi,colnames(RESP[,c(21:30+2)]))

# 2 factor model
# 確認的IRTモデルをおこなうためには，mirt.model関数で，モデルファイルを作っておく必要がある。
models <- mirt.model(
  "F1 = 1-25
   F2 = 1-25
   (F1*F2)=1-25
   COV=F1*F2
  "
)
mod2 <- mirt(ngaku16, model = models, itemtype = "2PL")

# 探索的IRTモデルの場合は，数値を指定するだけで良いそう。
mod2 <- mirt(ngaku16, model = 2, itemtype = "2PL", method = "EM", optimizer = "NR")

summary(mod2)
# これだとカテゴリカル因子分析の因子負荷量を計算してしまう。

# カテゴリカル因子分析のパラメタから，多次元IRTのパラメタへ変換する。
para <- coef(mod2,simplify=T)$items
plot(para[,1],para[,2])
MDIFF(mod2)
MDISC(mod2)
# MDIFFとMDISCを見る限り，一次元IRTモデルのパラメタとほぼ一緒のようだ

# practice for using "ltm" package
unidimTest(data = ngaku16)

system.time(
  res_2PL <- ltm(ngaku16 ~ z1)
)


res_3PL <- tpm(ngaku16)

# original function "MML_EM"
system.time(
  res <- MML_EM_cpp(ngaku16, fc=1)
)


# practice for using "irtoys" package
# bilog エンジンは仕えないもよう。

system.time(
  res_irtoys <- est(ngaku16, model = "2PL",engine = "ltm"  ,nqp = 31)
)

system.time(
  res_irtoys <- est(ngaku16, model = "2PL",engine = "icl"  ,nqp = 31)
)

# かなり早め

#practice for using "lazy.irt" package

cate <- rep(2,25)

system.time(
  res_lazy <- uIRT(ngaku16, type = rep("B2", 25))
)



# 2 parameter multidimensional logistic model

micc2pl <- function(theta,a,d,D = 1.702){
  # theta : multidimensional latent trait vector
  # a     : multidimensional discrimination parameter vector
  # b     : intercept in multidimensional space
  p <- 1/1+exp(-D*a%*%theta+d)
  return(p)
}

