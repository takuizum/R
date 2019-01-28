Likelihood & Bayes theorem
================
T.SHIBUYA
2018/10/12

# 尤度関数と事前分布，事後分布についての整理

## 尤度関数と最尤推定

``` r
dice_p <- function(w,t=3,s=10){
  
  w # サイコロで1がでる確率
  v <- 1-w # サイコロで1以外がでる確率
  t # あるひとつの目がでた回数
  s # 全部の試行の回数
  
  num <- prod(seq.int(s,s-t+1))
  den <- prod(seq.int(t:1))
  const <- num/den # 確率のための定数項，なくてもいい
  
  # 事後確率を求めるためにはサイコロすべての出目が分かる必要があるので，省略。
  res <- w^t * v^(s-t)* const # probability 
  return(res)
}
```

例えば，10回サイコロを振って，5回1がでたとする。その場合にサイコロの1の目がでる確率の密度関数は，先に定義したサイコロのパラメタの確率密度関数で定義可能であるが，実はこれは二項分布と同じである（1以外の確率と定義しているため，コインの裏表と同じ事をやっている）。

``` r
ggplot(data=data.frame(x=c(0:1)),aes(x=x))+
  stat_function(fun=dice_p,args = list(t=5,s=10))
```

![](20181012_Likelohood___Bayes_theorem_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

``` r
ggplot(data=data.frame((x=c(0:1))), aes(x=x))+
  stat_function(fun=dbinom, args=list(x=5,size=10))
```

![](20181012_Likelohood___Bayes_theorem_files/figure-gfm/unnamed-chunk-1-2.png)<!-- -->

この分布の形状は，正規化のための定数項を除いても変わらないため，以下の関数で定義してもよい。

``` r
dice_l <- function(w,t=3,s=10){
  
  w # サイコロで1がでる確率
  v <- 1-w # サイコロで1以外がでる確率
  
  res <- w^t * v^(s-t) # likelihood
  return(res)
}

ggplot(data=data.frame(x=c(0:1)),aes(x=x))+
  stat_function(fun=dice_l,args = list(t=5,s=10))
```

![](20181012_Likelohood___Bayes_theorem_files/figure-gfm/likelihood-1.png)<!-- -->

このように，定義される関数を尤度関数と呼ぶ。正規化のための定数項を取り除いたため，当然確率ではない。

さて，いまこのサイコロがどれくらい1の目がでやすいように偏っているかを知りたいとする。つまり知りたいのは1の目がでる確率（母数，パラメタ）である。これを知るためには実際にサイコロを投げてみて，そのデータにもとづいて推論をすればよい。例えば100回中35回，1の目がでた場合，

``` r
ggplot(data=data.frame(x=c(0:1)),aes(x=x))+
  stat_function(fun=dice_l,args = list(t=35,s=100))+
  geom_vline(xintercept = 35/100)
```

![](20181012_Likelohood___Bayes_theorem_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->
尤度関数のピーク（最尤推定値）は\(35/100\)のところに来ることが分かる。これは二項分布の最尤推定値がサンプル平均に一致するからである。指数分布族を使用すると期待値や分散，最尤推定値を，数値計算に依らず解析的に解くことができるので非常に便利。

このように尤度のピーク，最頻値をパラメタの推定値とする方法を最尤推定と呼ぶ。

## ベイズの定理との関係

ところで，たまたまサイコロを10回振ったときに1の目がでた場合に，そのサイコロのパラメタをサンプル平均の\(1/2\)とすることは，果たして妥当な推論だろうか。サイコロには6つの面があり，特にゆがみがなければひとつの目がでる確率は\(1/6\)であることは自明である。

しかし現実にはサイコロがちょっとゆがんでいるかもしれないし，たまたまそのような結果が得られる事もある。つまりサイコロの1の目がでるパラメタはある程度の幅を持った確率であると解釈する方がいいかもしれない。

ベイズ的な手法では，尤度以外に，サイコロの出目に関する事前の信念・情報を加え，データが得られた後のパラメタの分布を考える。このときの事前の信念・情報の確率分布を事前分布と呼び，データが得られた後のパラメタの確率分布を事後分布と呼ぶ。

いま，事前分布としてベータ分布をもちいる。この事前分布の形状を決定するためのパラメタをハイパーパラメタという。ハイパーパラメタは，1つ前の試行の結果や，これまでの先行研究の知見，常識などを考慮して選択される。今回は`shape1=2,
shape2=5`という設定を用いる。

``` r
dice_b <- function(w,t=3,s=10,a=2,b=5){
  
  w # サイコロで1がでる確率
  v <- 1-w # サイコロで1以外がでる確率
  
  res <- w^t * v^(s-t) * dbeta(w,a,b)# posterior ∝ likelihood * prior
  return(res)
}

ggplot(data=data.frame(x=c(0:1)),aes(x=x))+
  stat_function(fun=dbeta,args = list(shape1=2,shape2=5))
```

![](20181012_Likelohood___Bayes_theorem_files/figure-gfm/bayes-1.png)<!-- -->

この事前分布の情報を尤度関数に加えたときの事後分布の形状は，

``` r
ggplot(data=data.frame(x=c(0:1)),aes(x=x))+
  stat_function(fun=dice_b,args = list(t=5,s=10,a=2,b=5))
```

![](20181012_Likelohood___Bayes_theorem_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

のようになる。事前分布が正にゆがんだ分布であったため，事後分布の形状も左右対称の二項分布からすこし正にゆがんだ分布になる。

このとき，正規化のための定数項を入れていないことには注意されたい。