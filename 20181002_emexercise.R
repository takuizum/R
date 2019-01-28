# The standard error of estimate and information

# recovering item response matrix from table2.1 in text.
dat_in_ayala = matrix(ncol=5)
dat_in_ayala = numeric(0)

p_mat = matrix(nrow=2^5,ncol=5)
p = c(0,1)
j = 0
for(i in p){
  for(ii in p){
    for(iii in p){
      for(iiii in p){
        for(iiiii in p){
          j = j+1
          p_mat[j,] = c(iiiii,iiii,iii,ii,i)
        }
      }
    }
  }
}

p_mat = p_mat[p_mat %>% rowSums() %>% order(),]

freq = c(691,2280,242,235,158,184,1685,1053,134,462,92,65,571,79,87,41,1682,702,370,63,626,412,166,52,28,15,2095,1219,500,187,40,3385)
j = 0
for(n in freq){
  j = j + 1
  temp = matrix(p_mat[j,],nrow=n,ncol=5, byrow = T)
  dat_in_ayala = rbind(dat_in_ayala,temp)
}

devtools::install_github("takuizum/irtfun2")
library(irtfun2)

res_in_ayala = dat_in_ayala %>% irtfun2::estip(model="1PL", fc=1, EM_dist = 0, D=1)

theta_in_ayala = dat_in_ayala %>% irtfun2::estheta(param=res_in_ayala$para,est = "MLE", fc = 1, gc=0, method="NR")

res_in_ayala$para
res_in_ayala$SE
theta_in_ayala$res$MLE %>% unique()
theta_in_ayala$res$SE %>% unique()


# ICC

