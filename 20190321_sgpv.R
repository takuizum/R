# 20190321
############################################
# Second generation p-value
# https://github.com/weltybiostat/sgpv
############################################

devtools::install_github("weltybiostat/sgpv")
library(sgpv)

# Examples
## Simple example for three estimated log odds ratios but the same null interval
lb <- c(log(1.05), log(1.3), log(0.97))
ub <- c(log(1.8), log(1.8), log(1.02))
sgpv <- sgpvalue(est.lo = lb, est.hi = ub, null.lo = log(1/1.1), null.hi = log(1.1))
sgpv$p.delta
sgpv$delta.gap


?t.test
t_test_res <- t.test(1:10, y = c(7:20)) 
t_test_res$conf.int[1:2]

plot(extra ~ group, data = sleep)
with(sleep, t.test(extra[group == 1], extra[group == 2]))
t.test(extra ~ group, data = sleep)


# sgpv plot
data(leukstats)
plotsgpv(est.lo=leukstats$ci.lo, est.hi=leukstats$ci.hi,
         null.lo=-0.3, null.hi=0.3,
         set.order=order(leukstats$p.value),
         x.show=7000,
         plot.axis=c("TRUE","FALSE"),
         null.pt=0, outline.zone=TRUE,
         title.lab="Leukemia Example", y.lab="Fold Change (base 10)",
         x.lab="Classical p-value ranking",
         legend.on=TRUE)
axis(side=2,at=round(log(c(1/1000,1/100,1/10,1/2,1,2,10,100,1000),
                         base=10),2),labels=c("1/1000","1/100","1/10","1/2",1,2,10,100,1000),
     las=2)
