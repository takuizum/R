# ANOVA

seikaku <- factor(c(rep("外交的", 15), rep("内向的", 15)))
gakusyu <- factor(c(rep(c(rep("教師", 5),rep("教科書", 5), rep("タブレット", 5)), 2)))


score <- c(7,8,8,9,6,4,3,5,6,2,4,2,3,3,2,3,6,4,5,5,5,6,2,3,3,7,6,4,4,6)

# 分散分析を実行
res <- aov(score~seikaku*gakusyu)

summary(res)


# ANOVAくんで実行。
# ANOVAくんのテキストデータはこちらに
# http://riseki.php.xdomain.jp/index.php?plugin=attach&refer=ANOVA%E5%90%9B&openfile=anovakun_482.txt


anovakun(data.frame(seikaku,gakusyu, score), "ABs", 2,3)


b1 <- c(4.2,7.0,4.4)
b2 <- c(5.8,4.0,7.0)
lab <- c(1,2,3)
dat1 <- cbind(lab,b1)
dat2 <- cbind(lab,b2)

plot(dat1, type="b", ylim = c(0,8), xlim = c(0.5,3.5), xaxt="n",　xlab ="教材", ylab ="得点", pch = 16)
par(new = T)
plot(dat2, type="b", ylim = c(0,8), xlim = c(0.5,3.5), lty = 2, xaxt="n", xlab = "", ylab = "")
axis(1, c(1,2,3), c("教材A","教材B","教材C"))





