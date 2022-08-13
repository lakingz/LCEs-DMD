########################################################################
###    The Lyapunov Characteristic Exponents and their computation   ###
###                        By Ankai Liu                              ###
########################################################################
### Ref: The Lyapunov Characteristic Exponents and their computation ###
###                   By Charalampos Skokos                          ###
###                         page 44                                  ###
########################################################################

rm(list=ls())
library(deSolve)
library(reshape2)
library(ggplot2)
###########ODE decleration
DE_full <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS = (1-p)*u - b*S*I - u*S;
    dI = b*S*I - r*I - u*I;
    dS_l = (-b*I-u)*S_l - b*S*I_l;
    dI_l = b*I*S_l + (b*S-r-u)*I_l
    list(c(dS, dI, dS_l, dI_l))
  })
}

parameters <- c(p = 0.8, b = 200, u = 1/70, r = 20)
p <- parameters[1]
b <- parameters[2]
u <- parameters[3]
r <- parameters[4]
lambda1 <- -u; lambda1
lambda2 <- b*(1-p)-r-u; lambda2
namesofstates <- c('S','I','S_l','I_l')
initial <- c(S = 1, I = 1, S_l = 0, I_l = 0)
dim <- 2
####################
stepsize <- 0.01
begintime <- 0
endtime <- 500
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
#  ##Gram-schmidt (We do not use QR here since it will normlize the basis)
#  tangent <- state[,(dim+1):ncol(state)]
#  if (dim>1) {
#    for (i in c(2:dim)) {
#      for (k in c(1:(i-1))){
#        tangent[i,] <- tangent[i,]-tangent[k,]*sum(tangent[i,]*tangent[k,])
#      }
#    }  
#  }
#  state[,(dim+1):ncol(state)] <- tangent
  ##Apply svd to estimate dominating mode as well as growth rate
  
  
  
  ##
  for (i in c(1:dim)) {
    gamma[j,i+1] <- sqrt(sum(state[i,(dim+1):ncol(state)]*state[i,(dim+1):ncol(state)]))
    LE[j,i+1] <- sum(log(gamma[1:j,i+1]))/(j*itera_length)
    state[i,(dim+1):ncol(state)] <- state[i,(dim+1):ncol(state)]/gamma[j,i+1]
  }
}
##plot
LE_long <- melt(LE, id = 't')
ggplot(LE_long, aes(x = t,y = value, color = variable)) + geom_line()
tail(LE$`1`)
tail(LE$`2`)


###test bench

period_full <- seq(begintime+discard_time,endtime, by =stepsize)
LE_validation <- data.frame(t = period_full)
for (i in c(1:dim)) {
  initial_validation <- initial
  initial_validation[(dim+1):length(initial)] <- as.numeric(state[i,(dim+1):ncol(state)])
  out <- ode(y = initial_validation, times = period_full, func = DE_full, parms = parameters)
  ##Here we log every thing first to avoid 'inf'
  tangent_validation <- log(abs(out[,-c(1:(dim+1))]))
  lambda_tamp <- (2*tangent_validation[,c(1:(dim))]) %*% rep(1,dim)/out[,1]

  print(i)
  print(tangent_validation[,1])
  print(lambda_tamp)
  LE_validation <- cbind(LE_validation, lambda_tamp)  
}

colnames(LE_validation) <- c('t',1:dim)
head(LE_validation)
tail(LE_validation)
