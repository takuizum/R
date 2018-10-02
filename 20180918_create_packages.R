install.packages(c("devtools","usethis","roxygen2","testthat"))

library(devtools)
library(usethis)
library(roxygen2)
library(testthat)


create_package("irtfun")

use_mit_license("irtfun")

use_roxygen_md()

use_package_doc()

use_rcpp()


check()
