01\_introduction\_to\_irtfun2
================
T.SHIBUYA
2019/1/18

## `irtfun2` とは。

`irtfun2`は私が作成，維持している項目反応理論 (Item Response Theory, IRT)
の分析に特化したパッケージです。まだまだ機能は少ないですが，Rで使用可能なパッケージにはないいくつかの機能を備えています。

  - 多母集団モデルでのパラメタ推定ができる。
  - 項目パラメタの推定方法として，周辺ベイズ推定法や正則化周辺最尤推定法が使用できる。
  - 一般項目反応モデル (Genelal IRT model) の分析ができる。

などです。

　とはいえ，まずはソフトを走らせてみることが，理解への早道だと思いますので，早速`irtfun2`をインストールして，関数を使用してみましょう。

## `estip`関数

    starting httpd help server ... done

`Examples`のコードは仮想的に発生させた01のテストデータを`estip`で分析するものです。これを実行してみましょう。

``` r
res <- estip(x=sim_data_1, fc=2)
```

    Checking the model string vector 1 !
    The number of subject is 3000.
    The number of item is 30.
    The number of remove item is 0.
    Start calculating estimated item parameters.
    1 times -2 Marginal Log Likelihood is 77757.661287 
    2 times -2 Marginal Log Likelihood is 76738.402426 
    3 times -2 Marginal Log Likelihood is 76705.093466 
    4 times -2 Marginal Log Likelihood is 76703.139743 
    5 times -2 Marginal Log Likelihood is 76702.891910 
    6 times -2 Marginal Log Likelihood is 76702.849155 
    7 times -2 Marginal Log Likelihood is 76702.841027 
    8 times -2 Marginal Log Likelihood is 76702.839433 
    9 times -2 Marginal Log Likelihood is 76702.839117 
    10 times -2 Marginal Log Likelihood is 76702.839054 
    EM algorithm has been converged.
    Total iteration time is 10
    Start calculating estimated population distribution.
    1 times -2 Marginal Log Likelihood is 82778.891842 
    2 times -2 Marginal Log Likelihood is 76754.987282 
    3 times -2 Marginal Log Likelihood is 76695.641429 
    4 times -2 Marginal Log Likelihood is 76692.679855 
    5 times -2 Marginal Log Likelihood is 76692.181573 
    6 times -2 Marginal Log Likelihood is 76691.989642
