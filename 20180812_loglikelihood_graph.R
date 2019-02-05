# 尤度の考え方 in IRT

library(tidyverse);library(magick)

ptheta <- function(theta,a,b,c){
  c+(1-c)/(1+exp(-1.702*a*(theta-b)))
}

LL <- function(u,theta,a,b,c){
  sum(log(ptheta(theta,a,b,c))*u+log(1-ptheta(theta,a,b,c))*(1-u),na.rm = T)
}
LL_b <- function(u,theta,a,b,c,mu,sigma){
  sum(log(ptheta(theta,a,b,c))*u+log(1-ptheta(theta,a,b,c))*(1-u)-0.5*((theta-mu)/sigma)^2,na.rm = T)
}

a <- c(0.7,1.3,2.3,2.9)
b <- c(-2,-0.8,0,1.9)
c <- c(0.1,0.34,0.2,0.2)


x <- matrix(nrow=3^4,ncol=4)
pattern <- c(0,1,NA)
h <- 1
for(i in pattern){
  for(j in pattern){
    for(k in pattern){
      for(l in pattern){
        x[h,1] <- i
        x[h,2] <- j
        x[h,3] <- k
        x[h,4] <- l
        h <- h+1
      }
    }
  }
}

x <- x[order(rowSums(x,na.rm = T)),]
total <- rowSums(x,na.rm=T)


# 尤度関数
theta <- seq(-6,6,length.out = 301)
LOL <- matrix(nrow=301,ncol=3^4) %>% data.frame()

for(i in 1:3^4){
  LOL[,i] <- theta %>% matrix(ncol=1) %>% apply(1,LL,u=x[i,],a=a,b=b,c=c)
  colnames(LOL)[i] <- paste(x[i,1],x[i,2],x[i,3],x[i,4])
}
LOL$theta <- theta

LOL %>% tidyr::gather(key=pattern,value=lol,-theta) %>% 
  ggplot(aes(x=theta,y=lol,colour=pattern)) +
  geom_line()

# 尤度関数*事前分布
theta <- seq(-6,6,length.out = 301)
LOL_b <- matrix(nrow=301,ncol=3^4) %>% data.frame()

for(i in 1:3^4){
  LOL_b[,i] <- theta %>% matrix(ncol=1) %>% apply(1,LL_b,u=x[i,],a=a,b=b,c=c,mu=0,sigma=1)
  colnames(LOL_b)[i] <- paste(x[i,1],x[i,2],x[i,3],x[i,4])
}
LOL_b$theta <- theta

LOL_b %>% tidyr::gather(key=pattern,value=lol,-theta) %>% 
  ggplot(aes(x=theta,y=lol,colour=pattern)) +
  geom_line()


# 
mle <- estheta(cbind(c(1:81),x),param=cbind(a,b,c),est="MLE",ITEMc=2,Gc=0, method = "BFGS")
map <- estheta(cbind(c(1:81),x),param=cbind(a,b,c),est="MAP",ITEMc=2,Gc=0, method = "BFGS")
eap <- estheta(cbind(c(1:81),x),param=cbind(a,b,c),est="EAP",ITEMc=2,Gc=0)

optimise(LL,c(-10,10), maximum = T,u=x[37,],a=a,b=b,c=c)

optim(par=c(0),fn=LL,gr=function(u, theta,a,b,c){
  D <- 1.702
  p <- ptheta(theta,a,b,c)
  D*sum(a*(u - p)*(p-c)/(p*(1-c)), na.rm = T)
},method="BFGS",u=x[37,],control=list(fnscale=-1),a=a,b=b,c=c)

optim(par=c(0),fn=LL,method="SANN",u=x[37,],control=list(fnscale=-1),a=a,b=b,c=c)



img <- image_graph()

for(i in 1:3^4){
  cat("subject ",i,"\r")
  main <- sum(x[i,],na.rm = T)
  plot(x=LOL$theta,y=LOL[,i], ylim=c(-45,1),xlim=c(-6,6),ylab="lol",xlab=paste0("score is ",main,": pattern is ",colnames(LOL)[i]),type="l",lwd=2)
  par(new=T)
  plot(x=LOL_b$theta,y=LOL_b[,i], ylim=c(-45,1),xlim=c(-6,6),ylab="",xlab="",type="l",lty=2,lwd=2)
  mle_i <- mle$res$MLE[i]
  map_i <- map$res$MAP[i]
  eap_i <- eap$res$EAP[i]
  if(!is.na(mle_i)){
    par(new=T)
    plot(x=mle_i,y=LL(x[i,],mle_i,a,b,c), ylim=c(-45,1),xlim=c(-6,6),ylab="",xlab="",col="red", xaxt="n", yaxt="n")
    mtext(text=paste("MLE is",mle_i),line=0,col=2)
  } else {
    mtext(text="MLE is NA",line=0,col=2)
  }
  par(new=T)
  plot(x=map_i,y=LL_b(x[i,],map_i,a,b,c,mu=0,sigma=1), ylim=c(-45,1),xlim=c(-6,6),ylab="",xlab="",col=3, xaxt="n", yaxt="n")
  mtext(text=paste("MAP is",map_i),line=1,col=3)
  par(new=T)
  plot(x=eap_i,y=LL_b(x[i,],eap_i,a,b,c,mu=0,sigma=1), ylim=c(-45,1),xlim=c(-6,6),ylab="",xlab="",col=4, xaxt="n", yaxt="n")
  mtext(text=paste("EAP is",eap_i),line=2,col=4)
}

dev.off()
  
lol_anime <- image_animate(img, fps=5, loop=0)
image_write(lol_anime,path="lolanime.gif",format = "gif")

res_theta <- theta %>% data.frame()
res_theta$LL <- theta
for(p in 1:301){
  pp <- LL(x[37,],theta[p],a=a,b=b,c=c)
  res_theta[p,2] <- pp
}
