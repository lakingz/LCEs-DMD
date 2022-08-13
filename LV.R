########################################################################
###    The Lyapunov Characteristic Exponents and their computation   ###
###                        By Ankai Liu                              ###
########################################################################
### Ref: The Lyapunov Characteristic Exponents and their computation ###
###                   By Charalampos Skokos                          ###
###                         page 44                                  ###
########################################################################


rm(list=ls())
dev.off()
library(deSolve)
library(reshape2)
library(ggplot2)
###########ODE decleration


DE_full <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a * X + Y * Z
    dY <- b * (Y - Z)
    dZ <- -X * Y + c * Y - Z
    dX_l <- a * X_l + Y_l * Z + Y * Z_l
    dY_l <- b * (Y_l - Z_l)
    dZ_l <- -X_l * Y - X * Y_l + c * Y_l - Z_l
    list(c(dX, dY, dZ, dX_l, dY_l, dZ_l))
  })
}
parameters <- c(a = -8/3, b = -10, c = 28)
namesofstates <- c('X','Y','Z','X_l','Y_l','Z_l')
initial <- c(X = 1, Y = 1, Z = 1, X_l = 0, Y_l = 0, Z_l = 0)
dim <- 3
#out <- ode(y = initial, times = times, func = Lorenz_full, parms = parameters)
####################
stepsize <- 0.01
begintime <- 0
endtime <- 100
itera_length <- 0.1
itera_num <- endtime/itera_length
itera_period <- seq(begintime,itera_length, by = stepsize)
discard_time <- 100 ## we discard time interval up to t=100
discard_period <- seq(begintime,begintime+discard_time, by =stepsize)
out <- ode(y = initial, times = discard_period, func = DE_full, parms = parameters)
initial <- out[nrow(out),-1]
#####################memory allocation
state <- data.frame(t(initial))
gamma <- data.frame(t= c(1:itera_num)*itera_length)
LE <- data.frame(t= c(1:itera_num)*itera_length)
for (i in c(1:dim)) {
  state <- rbind(state,t(initial))  
  state[i+1,dim+i] <- state[i+1,dim+i] + 1
  
  gamma <- cbind(gamma, NA)
  colnames(gamma)[i+1] <- i
  LE <- cbind(LE, NA)
  colnames(LE)[i+1] <- i
}
state <- state[-1,] 
rownames(state) <- c()
##main iteration
for (j in c(1:itera_num)) {
  orbit <- state[1,1:dim]
  for (i in c(1:dim)) {
    
    state[i,1:dim] <- orbit  #we force to stay on the same orbit
    out <- ode(y = setNames(as.numeric(state[i,]),namesofstates), 
               times = begintime+itera_length*(j-1)+itera_period, func = DE_full, parms = parameters)
    state[i,] <- out[nrow(out),-1]
  }
  ##Gram-schmidt (We do not use QR here since it will normlize the basis)
  tangent <- state[,(dim+1):ncol(state)]
  if (dim>1) {
    for (i in c(2:dim)) {
      for (k in c(1:(i-1))){
        tangent[i,] <- tangent[i,]-tangent[k,]*sum(tangent[i,]*tangent[k,])
      }
    }  
  }
  state[,(dim+1):ncol(state)] <- tangent
  ##
  for (i in c(1:dim)) {
    gamma[j,i+1] <- sqrt(sum(state[i,(dim+1):ncol(state)]*state[i,(dim+1):ncol(state)]))
    LE[j,i+1] <- sum(log(gamma[1:j,i+1]))/(j*itera_length)
    state[i,(dim+1):ncol(state)] <- state[i,(dim+1):ncol(state)]/gamma[j,i+1]
  }
}
##plot
tail(LE)
LE_long <- melt(LE, id = 't')
ggplot(LE_long, aes(x = t,y = value, color = variable)) + geom_line()
