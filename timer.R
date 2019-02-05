# stop whatch

#time <- 0
#repeat {
#  if(time == 0) START <- Sys.time()
#  t0 <- Sys.time()
#  
#  DIFF <- 0
# while(DIFF == 0){
#    
#    t1 <- Sys.time()
#    
#    DIFF <- t1 - t0
#    
#    if(DIFF < 0.01) DIFF <- 0
#    
#  }
  
#  time <- t1 - START
#  
#  time <- round(as.numeric(time, units = "secs"), digits = 3)
#  cat(time, "sec.\r")
#}









# functionalize







timer <- function(StopTime, MESSAGE = "STOP!!", disp = T, units = "secs", 
                  cex = 5, # だいたい10位がベスト 
                  font = 1, 
                  col = 1, 
                  family="",  # Macの場合は，Helveticaを選択すればOK，ただし日本語は不可。
                  os = "win"
                  ){
  
  time <- 0
  if(disp == T){
    
    repeat {
      if(time == 0) START <- Sys.time()
      t0 <- Sys.time()
      
      DIFF <- 0
      while(DIFF == 0){
        
        t1 <- Sys.time()
        
        DIFF <- t1 - t0
        
        if(DIFF < 0.01) DIFF <- 0
        
      }
      
      time <- t1 - START
      
      if(time >= StopTime) break
      
      time <- round(as.numeric(time, units = units), digits = 3)
      cat(time, units,".\r")
    }
    
    cat("\n",MESSAGE)
    
    if(os =="win")    win.graph()
    if(os == "mac") quartz()
    
    plot(NA, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n")
    par(cex =cex, col = col, font = font, family = family)
    text(0.5,0.5,MESSAGE)
    
  }
 
  
}
