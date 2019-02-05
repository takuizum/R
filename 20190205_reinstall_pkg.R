# installation at the new env

# essentials
need_packages <- c("Rcpp", "tidyverse", "rstan")
install.packages(need_packages)

# always use
need_pkg2 <- c("pracma", "psych", "irtoys", "ltm", "equate")
install.packages(need_pkg2)
# install.packages("sm")

# useful
need_pkg3 <- c("magrittr", "devtools", "shinystan")
install.packages(need_pkg3)

# git version
devtools::install_github("tidyverse/ggplot2")

# irt tools
need_pkg4 <- c("plink", "mirt", )


library(irtoys)

?est
p.1pl <- est(Scored, model="1PL", engine="icl")

# irtfun2 sourse ver
devtools::install_github("takuizum/irtfun2", type = "sourse", dependencies = TRUE)
