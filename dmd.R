######################################DMD
data <- as.matrix(out[,-1])
A_tilda <- d_inv%*%t(u)%*%data_out%*%v
A_tilda
ev <- eigen(t(A_tilda))
eig <- t(ev$vectors)
lambda <- ev$values
lambda
phi[,1]%*%d_inv%*%t(u)%*%data_out
############################
data <- input$VORTALL
data_in <- data[,-1]
data_out <- data[,-ncol(data)]

svd_decomp <- svd(data_in)
u <- svd_decomp$u
v <- svd_decomp$v
d <- diag(svd_decomp$d)
d_inv <- diag(1/svd_decomp$d)

A_tilda <- t(u)%*%data_out%*%v%*%d_inv
A_tilda
ev <- eigen(A_tilda)
eig <- ev$vectors
lambda <- ev$values
lambda
phi <- data_out%*%v%*%d_inv%*%eig
phi
dev.off()

plot(phi[3,],type = 'p')
####u multiply from the left while v multiply from the right.
#t(u)%*%u
#v%*%t(v)



##########parameter for breeding
bred_cycle <- 10
initial_bred_radius <- 0.01
sample_points <- 10
initial_per_cycle <- matrix(rnorm(sample_points*length(state),
                                  mean = state,sd = initial_bred_radius),
                            ncol = length(state))
initial_per_cycle <- rbind(state,initial_per_cycle)


for (i in 1:bred_cycle) {
  out <- ode(func = Parasite, y = c(P = 0.5, H = 0.5), times = 0:50, parms = ks,
             method = "iteration")
  ?apply
}
out <- ode(func = Parasite, y = c(P = 0.5, H = 0.5), times = 0:50, parms = ks,
           method = "iteration")

out2<- ode(func = Parasite, y = c(P = 0.5, H = 0.5), times = 0:50, parms = 25,
           method = "iteration")

out3<- ode(func = Parasite, y = c(P = 0.5, H = 0.5), times = 0:50, parms = 35,
           method = "iteration")

## Plot all 3 scenarios in one figure
plot(out, out2, out3, lty = 1, lwd = 2)

