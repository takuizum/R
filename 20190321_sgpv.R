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
