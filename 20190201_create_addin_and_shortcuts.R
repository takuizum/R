# 
devtools::install_github("rstudio/addinexamples", type = "binary")

# Ctrl + Shift + I  = %in%
#  %in% 

#
library(rstudioapi)

insertPipeinAddin <- function(){
  rstudioapi::insertText(" %<>% ")
}
  
  
dat <- c(1,2,3,4)
dat %>% data.frame()


dat %<>% data.frame()
dat
