# observed score equating

# equi percentile equating


# gen data
dat <- equate::ACTmath

X <- as.freqtab(ACTmath[, 1:2])
Y <- as.freqtab(ACTmath[, c(1, 3)])


# equate package

res_equate <- equate::equate(x = X, y = Y, type = "equipercentile")
res_equate$concordance


# irtfun2

X <- dat$scale %>% rep(dat$xcount)
Y <- dat$scale %>% rep(dat$ycount)

irtfun2::epe(x = X, y = Y)
