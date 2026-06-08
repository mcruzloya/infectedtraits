# Vector competence thermal response analysis
# Fits semi-mechanistic models for infection and dissemination probability 
# using JAGS.
# Run from the repository root directory.

set.seed(42) # Set seed for reproducibility.

# We will need these packages. Please make sure you install them beforehand!
library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

# FlexTPC model for viral replication rate.
flexTPC <- function(T, Tmin, Tmax, rmax, alpha, beta) {
  s <- alpha * (1 - alpha) / beta^2
  result <- rep(0, length(T))
  Tidx = (T > Tmin) & (T < Tmax)
  result[Tidx] <- rmax * exp(s * (alpha * log( (T[Tidx] - Tmin) / alpha) 
                                  + (1 - alpha) * log( (Tmax - T[Tidx]) / (1 - alpha))
                                  - log(Tmax - Tmin)) ) 
  return(result)
}

# Hill function for immunity
Hill <- function(T, rmin, rmax, T50, n) {
  a = pmin(exp(n * log(T)), 10^100)
  b = pmin(exp(n * log(T50)), 10^100)
  return(pmax(rmin + (rmax - rmin) * a / (a + b), 0.0))
}

# prop. infected model
bbc <- function(T, im_rmin, im_rmax, im_T50, im_n,
                vr_Tmin, vr_Tmax, vr_rmax, vr_alpha, vr_beta,
                bbc_n) {
  im <- Hill(T, im_rmin, im_rmax, im_T50, im_n)
  vr <- flexTPC(T, vr_Tmin, vr_Tmax, vr_rmax, vr_alpha, vr_beta)
  
  v <- rep(0, length(T))
  idx <- (vr > 0.0) & (vr > im)
  
  v[idx] <- 1.0 - im[idx] / vr[idx]
  
  bbc <- exp(bbc_n * log(v))
  return(bbc)
}

data <- read.csv("./data/vector_competence.csv")
data$p.infected.sderr <- sqrt((data$p.infected * (1 - data$p.infected)) / data$total)
data$p.disseminated.sderr <- sqrt((data$p.disseminated * (1 - data$p.disseminated)) / data$total)

data.WN02 <- subset(data, (data$genotype == "WN02") & (data$cycle == FALSE))
data.NY10 <- subset(data, (data$genotype == "NY10") & (data$cycle == FALSE))

data.WN02

errbar.length=0.06
errbar.width=2

temps <- seq(0, 40, 0.1)

par(mfrow=c(2,2))

plot(data.WN02$temperature, data.WN02$p.infected, main="WN02",
     xlab="Temperature [°C]", ylab="b (prop. infected)", xlim=c(0, 40),
     ylim=c(0,1), pch=20)
arrows(data.WN02$temperature, data.WN02$p.infected, y1=data.WN02$p.infected+data.WN02$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.WN02$temperature, data.WN02$p.infected, y1=data.WN02$p.infected-data.WN02$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)


plot(data.NY10$temperature, data.NY10$p.infected, main="NY10",
     xlab="Temperature [°C]", ylab="b (prop. infected)", xlim=c(0, 40),
     ylim=c(0,1), pch=20)
arrows(data.NY10$temperature, data.NY10$p.infected, y1=data.NY10$p.infected+data.NY10$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.NY10$temperature, data.NY10$p.infected, y1=data.NY10$p.infected-data.NY10$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

plot(data.WN02$temperature, data.WN02$p.disseminated, main="WN02",
     xlab="Temperature [°C]", ylab="bc (prop. disseminated)", xlim=c(0, 40),
     ylim=c(0,1), pch=20)
arrows(data.WN02$temperature, data.WN02$p.disseminated, y1=data.WN02$p.disseminated+data.WN02$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.WN02$temperature, data.WN02$p.disseminated, y1=data.WN02$p.disseminated-data.WN02$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

plot(data.NY10$temperature, data.NY10$p.disseminated, main="NY10",
     xlab="Temperature [°C]", ylab="bc (prop. disseminated)", xlim=c(0, 40),
     ylim=c(0,1), pch=20)
arrows(data.NY10$temperature, data.NY10$p.disseminated, y1=data.NY10$p.disseminated+data.NY10$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.NY10$temperature, data.NY10$p.disseminated, y1=data.NY10$p.disseminated-data.NY10$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

sink("bbc_mech.txt")
cat("
    model{
    ## Priors for mosquito viral immunity model.
    im_T50 ~ dnorm(20, 1 / 10^2)
    im_logn ~ dnorm(0, 1 / 0.5^2)
    im_n <- 10^im_logn
    im_rmin ~ dunif(0, 1) # Probably small but nonzero but setting uniform prior.
    im_rmax <- 1.0 # Normalized to one for maximum immunity.
    
    ## Priors for WN02 viral replication rate.
    vr_WN02_Tmin ~ dnorm(5, 1 / 2.5^2)
    vr_WN02_Tmax ~ dnorm(35, 1 / 5^2)
    vr_WN02_rmax ~ dexp(1 / 1)
    vr_WN02_alpha <- 0.8
    vr_WN02_beta <- 0.2
    
    # Transformed parameters.
    vr_WN02_s <- vr_WN02_alpha * (1 - vr_WN02_alpha) / vr_WN02_beta^2
    vr_WN02_Topt <- vr_WN02_alpha * vr_WN02_Tmax + (1 - vr_WN02_alpha) * vr_WN02_Tmin
    
    ## Priors for NY10 viral replication rate.
    vr_NY10_Tmin ~ dnorm(5, 1 / 2.5^2)
    vr_NY10_Tmax ~ dnorm(35, 1 / 5^2)
    vr_NY10_rmax ~ dexp(1 / 1)
    vr_NY10_alpha <- 0.8
    vr_NY10_beta <- 0.2
    
    # Transformed parameters
    vr_NY10_s <- vr_NY10_alpha * (1 - vr_NY10_alpha) / vr_NY10_beta^2
    vr_NY10_Topt <- vr_NY10_alpha * vr_NY10_Tmax + (1 - vr_NY10_alpha) * vr_NY10_Tmin

    ## Priors for proportion infected model
    b_WN02_v50 ~ dbeta(0.2*30, 0.8*30)
    b_WN02_n <- -log(2) / log(b_WN02_v50)
    
    b_NY10_v50 ~ dbeta(0.2*30, 0.8*30)
    b_NY10_n <- -log(2) / log(b_NY10_v50)
    
    ## Priors for proportion disseminated model
    bc_WN02_v50 ~ dbeta(0.3*30, 0.7*30)T(b_WN02_v50, )
    bc_WN02_n <- -log(2) / log(bc_WN02_v50)

    bc_NY10_v50 ~ dbeta(0.3*30, 0.7*30)T(b_NY10_v50, )
    bc_NY10_n <- -log(2) / log(bc_NY10_v50)

    ## Likelihood
    for(i in 1:N.temps){
      # Immunity model (Hill equation).
      im[i] <- im_rmin + (im_rmax - im_rmin) * (temp[i]^im_n /
                                              (im_T50^im_n + temp[i]^im_n) )
    
      # Viral replication rate models (restritcted flexTPC equation).
      vr_WN02[i] <- (vr_WN02_Tmax > temp[i]) * (vr_WN02_Tmin < temp[i]) * vr_WN02_rmax * exp(vr_WN02_s * (
                      - log(vr_WN02_Tmax - vr_WN02_Tmin) 
                      + vr_WN02_alpha * log( max((temp[i] - vr_WN02_Tmin) / vr_WN02_alpha, 10^-15)) 
                      + (1 - vr_WN02_alpha) * log( max((vr_WN02_Tmax - temp[i]) / (1 - vr_WN02_alpha), 10^-15))))
      
      vr_NY10[i] <- (vr_NY10_Tmax > temp[i]) * (vr_NY10_Tmin < temp[i]) * vr_NY10_rmax * exp(vr_NY10_s * (
                      - log(vr_NY10_Tmax - vr_NY10_Tmin) 
                      + vr_NY10_alpha * log( max((temp[i] - vr_NY10_Tmin) / vr_NY10_alpha, 10^-15)) 
                      + (1 - vr_NY10_alpha) * log( max((vr_NY10_Tmax - temp[i]) / (1 - vr_NY10_alpha), 10^-15))))
      
      # Steady state viral load.
      v_WN02[i] <- ifelse(vr_WN02[i] == 0.0, 0.0, 
                   ifelse(vr_WN02[i] <= im[i], 0.0, 1.0 - im[i] / vr_WN02[i])
                   )
                   
      v_NY10[i] <- ifelse(vr_NY10[i] == 0.0, 0.0, 
                   ifelse(vr_NY10[i] <= im[i], 0.0, 1.0 - im[i] / vr_NY10[i])
                   ) 
      
      # Model for infection and dissemination proportions.
      b_WN02[i] <- max(v_WN02[i]^b_WN02_n, 10^-10)
      bc_WN02[i] <- max(v_WN02[i]^bc_WN02_n, 10^-10)
      
      b_NY10[i] <- max(v_NY10[i]^b_NY10_n, 10^-10)
      bc_NY10[i] <- max(v_NY10[i]^bc_NY10_n, 10^-10)
      
      n_b_WN02[i] ~ dbin(b_WN02[i], N_WN02[i])
      n_bc_WN02[i] ~ dbin(bc_WN02[i], N_WN02[i])
      
      n_b_NY10[i] ~ dbin(b_NY10[i], N_NY10[i])
      n_bc_NY10[i] ~ dbin(bc_NY10[i], N_NY10[i])
    }
    } # close model
    ",fill=TRUE)
sink()





##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  im_T50 = runif(1, min=10, max=20),
  im_logn = runif(1, min=0, max=1),
  im_rmin = runif(1, min=0, max=0.2),
  
  vr_WN02_Tmin = runif(1, min=0, max=10),
  vr_WN02_Tmax = runif(1, min=35, max=40),
  vr_WN02_rmax = runif(1, min=1, max=3),
  
  vr_NY10_Tmin = runif(1, min=0, max=10),
  vr_NY10_Tmax = runif(1, min=35, max=40),
  vr_NY10_rmax = runif(1, min=1, max=3),
  
  b_WN02_v50 = runif(1, min=0.1, max=0.3),
  b_NY10_v50 = runif(1, min=0.1, max=0.3)
  )}

##### Parameters to save in chain.
parameters <- c("im_T50", "im_logn", "im_n", "im_rmin", 
                "vr_WN02_Tmin", "vr_WN02_Tmax", "vr_WN02_rmax", "vr_WN02_Topt", "vr_WN02_alpha", "vr_WN02_beta",
                "vr_NY10_Tmin", "vr_NY10_Tmax", "vr_NY10_rmax", "vr_NY10_Topt", "vr_NY10_alpha", "vr_NY10_beta",
                "b_WN02_v50", "b_WN02_n",
                "bc_WN02_v50", "bc_WN02_n",
                "b_NY10_v50", "b_NY10_n",
                "bc_NY10_v50", "bc_NY10_n")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 2800000 # number of iterations in each chain
nb <- 800000 # number of 'burn in' iterations to discard
nt <- 64 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains

data.WN02
##### Organize Data for JAGS
jag.data <- list(temp=data.WN02$temperature, 
                 N.temps=length(data.WN02$total),
                 
                 N_WN02=data.WN02$total, 
                 n_b_WN02=data.WN02$n.infected,
                 n_bc_WN02=data.WN02$n.disseminated,
                 
                 N_NY10=data.NY10$total, 
                 n_b_NY10=data.NY10$n.infected,
                 n_bc_NY10=data.NY10$n.disseminated
                 )

##### Run JAGS
#mech.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                    model.file="bbc_mech.txt", n.thin=nt, n.chains=nc, 
#                    n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
mech.out <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                        model.file="bbc_mech.txt", 
                                        n.thin=nt, n.chains=nc, 
                                        n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                        jags.seed = 42))

mech.out
mcmcplot(mech.out)

par(mfrow=c(2,2))


chains.WN02.inf <- MCMCchains(mech.out, params=c("im_T50", "im_n", "im_rmin", 
                                                 "vr_WN02_Tmin", "vr_WN02_Tmax", "vr_WN02_rmax", "vr_WN02_alpha", "vr_WN02_beta",
                                                 "b_WN02_v50", "b_WN02_n"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.WN02.inf, 1, 
                function(x) bbc(temps, x[["im_rmin"]], 1.0, x[["im_T50"]], x[["im_n"]], 
                                       x[["vr_WN02_Tmin"]], x[["vr_WN02_Tmax"]], x[["vr_WN02_rmax"]], x[["vr_WN02_alpha"]], x[["vr_WN02_beta"]],
                                      x[["b_WN02_n"]]))

# Find mean curve and credible intervals.
meancurve.WN02.inf <- apply(curves, 1, mean)
CI.WN02.inf <- apply(curves, 1, quantile, c(0.025, 0.975))



plot(data.WN02$temperature, data.WN02$p.infected, main="WN02",
     xlab="Temperature [°C]", ylab="b (prop. infected)", xlim=c(0, 40),
     ylim=c(0,1), pch=20)
arrows(data.WN02$temperature, data.WN02$p.infected, y1=data.WN02$p.infected+data.WN02$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.WN02$temperature, data.WN02$p.infected, y1=data.WN02$p.infected-data.WN02$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.inf, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.WN02.inf[1,],
                                rev(CI.WN02.inf[2,])), 
        col=alpha("steelblue", 0.3), lty=0)


chains.NY10.inf <- MCMCchains(mech.out, params=c("im_T50", "im_n", "im_rmin", 
                                                 "vr_NY10_Tmin", "vr_NY10_Tmax", "vr_NY10_rmax", "vr_NY10_alpha", "vr_NY10_beta",
                                                 "b_NY10_n"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.NY10.inf, 1, 
                function(x) bbc(temps, x[["im_rmin"]], 1.0, x[["im_T50"]], x[["im_n"]], 
                                x[["vr_NY10_Tmin"]], x[["vr_NY10_Tmax"]], x[["vr_NY10_rmax"]], x[["vr_NY10_alpha"]], x[["vr_NY10_beta"]],
                                x[["b_NY10_n"]]))

# Find mean curve and credible intervals.
meancurve.NY10.inf <- apply(curves, 1, mean)
CI.NY10.inf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.NY10$temperature, data.NY10$p.infected, main="NY10",
     xlab="Temperature [°C]", ylab="b (prop. infected)", xlim=c(0, 40),
     ylim=c(0,1), pch=20)
arrows(data.NY10$temperature, data.NY10$p.infected, y1=data.NY10$p.infected+data.NY10$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.NY10$temperature, data.NY10$p.infected, y1=data.NY10$p.infected-data.NY10$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.inf, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.NY10.inf[1,],
                                rev(CI.NY10.inf[2,])), 
        col=alpha("steelblue", 0.3), lty=0)


chains.WN02.dis <- MCMCchains(mech.out, params=c("im_T50", "im_n", "im_rmin", 
                                                 "vr_WN02_Tmin", "vr_WN02_Tmax", "vr_WN02_rmax", "vr_WN02_alpha", "vr_WN02_beta",
                                                 "bc_WN02_n"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.WN02.dis, 1, 
                function(x) bbc(temps, x[["im_rmin"]], 1.0, x[["im_T50"]], x[["im_n"]], 
                                x[["vr_WN02_Tmin"]], x[["vr_WN02_Tmax"]], x[["vr_WN02_rmax"]], x[["vr_WN02_alpha"]], x[["vr_WN02_beta"]],
                                x[["bc_WN02_n"]]))

# Find mean curve and credible intervals.
meancurve.WN02.dis <- apply(curves, 1, mean)
CI.WN02.dis <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(data.WN02$temperature, data.WN02$p.disseminated, main="WN02",
     xlab="Temperature [°C]", ylab="bc (prop. disseminated)", xlim=c(0, 40),
     ylim=c(0, 1), pch=20)
arrows(data.WN02$temperature, data.WN02$p.disseminated, y1=data.WN02$p.disseminated +data.WN02$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.WN02$temperature, data.WN02$p.disseminated, y1=data.WN02$p.disseminated-data.WN02$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.dis, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.WN02.dis[1,],
                                rev(CI.WN02.dis[2,])), 
        col=alpha("steelblue", 0.3), lty=0)

chains.NY10.dis <- MCMCchains(mech.out, params=c("im_T50", "im_n", "im_rmin", 
                                                 "vr_NY10_Tmin", "vr_NY10_Tmax", "vr_NY10_rmax", "vr_NY10_alpha", "vr_NY10_beta",
                                                 "bc_NY10_n"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.NY10.dis, 1, 
                function(x) bbc(temps, x[["im_rmin"]], 1.0, x[["im_T50"]], x[["im_n"]], 
                                x[["vr_NY10_Tmin"]], x[["vr_NY10_Tmax"]], x[["vr_NY10_rmax"]],  x[["vr_NY10_alpha"]], x[["vr_NY10_beta"]],
                                x[["bc_NY10_n"]]))

# Find mean curve and credible intervals.
meancurve.NY10.dis <- apply(curves, 1, mean)
CI.NY10.dis <- apply(curves, 1, quantile, c(0.025, 0.975))


plot(data.NY10$temperature, data.NY10$p.disseminated, main="NY10",
     xlab="Temperature [°C]", ylab="bc (prop. disseminated)", xlim=c(0, 40),
     ylim=c(0,1), pch=20)
arrows(data.NY10$temperature, data.NY10$p.disseminated, y1=data.NY10$p.disseminated+data.NY10$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.NY10$temperature, data.NY10$p.disseminated, y1=data.NY10$p.disseminated-data.NY10$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)


# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.dis, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.NY10.dis[1,],
                                rev(CI.NY10.dis[2,])), 
        col=alpha("steelblue", 0.3), lty=0)


par(mfrow=c(1, 1))
temps <- seq(0, 40, 0.1) 

curves.im <- apply(chains.NY10.dis, 1, 
                   function(x) Hill(temps, x[["im_rmin"]], 1.0, x[["im_T50"]], x[["im_n"]]))
meancurve.im <- apply(curves.im, 1, mean)
CI.im <- apply(curves.im, 1, quantile, c(0.025, 0.975))


plot(temps, meancurve.im, col="gold", xlab="Temperature [°C]", 
     ylab="rate", xlim=c(0, 40),
     ylim=c(0, 5), type='l', main="Viral replication and mosquito immunity")
polygon(c(temps, rev(temps)), c(CI.im[1,],
                                rev(CI.im[2,])), 
        col=alpha("gold", 0.3), lty=0)

curves.vr.WN02 <- apply(chains.WN02.dis, 1, 
                   function(x) flexTPC(temps, x[["vr_WN02_Tmin"]], x[["vr_WN02_Tmax"]], x[["vr_WN02_rmax"]], x[["vr_WN02_alpha"]], x[["vr_WN02_beta"]]))
meancurve.vr.WN02 <- apply(curves.vr.WN02, 1, mean)
CI.vr.WN02 <- apply(curves.vr.WN02, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.vr.WN02, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.vr.WN02[1,],
                                rev(CI.vr.WN02[2,])), 
        col=alpha("darkgreen", 0.3), lty=0)

curves.vr.NY10 <- apply(chains.NY10.dis, 1, 
                        function(x) flexTPC(temps, x[["vr_NY10_Tmin"]], x[["vr_NY10_Tmax"]], x[["vr_NY10_rmax"]], x[["vr_NY10_alpha"]], x[["vr_NY10_beta"]]))
meancurve.vr.NY10 <- apply(curves.vr.NY10, 1, mean)
CI.vr.NY10 <- apply(curves.vr.NY10, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.vr.NY10, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.vr.NY10[1,],
                                rev(CI.vr.NY10[2,])), 
        col=alpha("darkblue", 0.3), lty=0)
legend(5, 4, c("immunity", "replication WN02", "replication NY10"), col=c("gold", "darkgreen", "darkblue"),
       lty=c(1, 1, 1), cex=0.8)


plot(temps, meancurve.im, col="gold", xlab="Temperature [°C]", 
     ylab="rate", xlim=c(5, 15),
     ylim=c(0, 1), type='l', main="Viral replication and mosquito immunity")
polygon(c(temps, rev(temps)), c(CI.im[1,],
                                rev(CI.im[2,])), 
        col=alpha("gold", 0.3), lty=0)

curves.vr.WN02 <- apply(chains.WN02.dis, 1, 
                        function(x) flexTPC(temps, x[["vr_WN02_Tmin"]], x[["vr_WN02_Tmax"]], x[["vr_WN02_rmax"]], x[["vr_WN02_alpha"]], x[["vr_WN02_beta"]]))
meancurve.vr.WN02 <- apply(curves.vr.WN02, 1, mean)
CI.vr.WN02 <- apply(curves.vr.WN02, 1, quantile, c(0.025, 0.975))

#lines(temps, meancurve.vr.WN02, col="darkgreen")
#polygon(c(temps, rev(temps)), c(CI.vr.WN02[1,],
#                                rev(CI.vr.WN02[2,])), 
#        col=alpha("darkgreen", 0.3), lty=0)

curves.vr.NY10 <- apply(chains.NY10.dis, 1, 
                        function(x) flexTPC(temps, x[["vr_NY10_Tmin"]], x[["vr_NY10_Tmax"]], x[["vr_NY10_rmax"]], x[["vr_NY10_alpha"]], x[["vr_NY10_beta"]]))
meancurve.vr.NY10 <- apply(curves.vr.NY10, 1, mean)
CI.vr.NY10 <- apply(curves.vr.NY10, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.vr.NY10, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.vr.NY10[1,],
                                rev(CI.vr.NY10[2,])), 
        col=alpha("darkblue", 0.3), lty=0)
legend(6, 0.8, c("immunity", "replication WN02", "replication NY10"), col=c("gold", "darkgreen", "darkblue"),
       lty=c(1, 1, 1), cex=0.8)


par(mfrow=c(1,2))
vir <- seq(0, 1, 0.02)
chains.WN02.bbc <- MCMCchains(mech.out, params=c("b_WN02_n",
                                                 "bc_WN02_n"))
b.curves <- apply(chains.WN02.bbc, 1, 
                  function(x)  vir^x[["b_WN02_n"]])
bc.curves <- apply(chains.WN02.bbc, 1, 
                  function(x)  vir^x[["bc_WN02_n"]])

meancurve <- apply(b.curves, 1, mean)
CI <- apply(b.curves, 1, quantile, c(0.025, 0.975))

plot(vir, meancurve, col="darkgreen", xlab="viral load", 
     ylab="probability", xlim=c(0, 1),
     ylim=c(-0.1, 1.1), type='l', main="WN02")
polygon(c(vir, rev(vir)), c(CI[1,],
                                rev(CI[2,])), 
        col=alpha("darkgreen", 0.3), lty=0)

meancurve <- apply(bc.curves, 1, mean)
CI <- apply(bc.curves, 1, quantile, c(0.025, 0.975))

lines(vir, meancurve, col="violet")
polygon(c(vir, rev(vir)), c(CI[1,],
                                rev(CI[2,])), 
        col=alpha("violet", 0.3), lty=0)
legend(0.35, 0.2, c("infection", "dissemination"), col=c("darkgreen", "violet"),
       lty=c(1, 1), cex=0.8)

#rmin, rmax, T50, n
chains.NY10.bbc <- MCMCchains(mech.out, params=c("b_NY10_n",
                                                 "bc_NY10_n"))


b.curves <- apply(chains.NY10.bbc, 1, 
                  function(x)  vir^x[["b_NY10_n"]])
bc.curves <- apply(chains.NY10.bbc, 1, 
                   function(x) vir^x[["bc_NY10_n"]])

meancurve <- apply(b.curves, 1, mean)
CI <- apply(b.curves, 1, quantile, c(0.025, 0.975))

plot(vir, meancurve, col="darkblue", xlab="viral load", 
     ylab="probability", xlim=c(0, 1),
     ylim=c(-0.1, 1.1), type='l', main="NY10")
polygon(c(vir, rev(vir)), c(CI[1,],
                            rev(CI[2,])), 
        col=alpha("darkblue", 0.3), lty=0)

meancurve <- apply(bc.curves, 1, mean)
CI <- apply(bc.curves, 1, quantile, c(0.025, 0.975))

lines(vir, meancurve, col="violet")
polygon(c(vir, rev(vir)), c(CI[1,],
                            rev(CI[2,])), 
        col=alpha("violet", 0.3), lty=0)
legend(0.35, 0.2, c("infection", "dissemination"), col=c("darkblue", "violet"),
       lty=c(1, 1), cex=0.8)
saveRDS(mech.out, "bbc_mech.RDS")
