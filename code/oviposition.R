# Oviposition thermal response analysis
# Fits flexTPC models to mosquito oviposition data using JAGS.
# Run from the repository root directory.

set.seed(42)

library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')


# FlexTPC model for thermal performance curves.
flexTPC <- function(T, Tmin, Tmax, rmax, alpha, beta) {
  s <- alpha * (1 - alpha) / beta^2
  result <- rep(0, length(T))
  Tidx = (T > Tmin) & (T < Tmax)
  result[Tidx] <- rmax * exp(s * (alpha * log( (T[Tidx] - Tmin) / alpha) 
                                  + (1 - alpha) * log( (Tmax - T[Tidx]) / (1 - alpha))
                                  - log(Tmax - Tmin)) ) 
  return(result)
}

data <- read.csv("./data/ovipositing.csv")
data
data$se <- sqrt(data$p.ovi * (1 - data$p.ovi) / data$total)


data.pO.unexp <- subset(data, (data$strain == "none") & (data$cycle == FALSE))
data.pO.WN02.uninf <- subset(data, (data$strain == "WN02") & (data$cycle == FALSE)
                             & (data$infected == FALSE))
data.pO.NY10.uninf <- subset(data, (data$strain == "NY10") & (data$cycle == FALSE)
                             & (data$infected == FALSE))
data.pO.WN02.inf <- subset(data, (data$strain == "WN02") & (data$cycle == FALSE)
                             & (data$infected == TRUE))
data.pO.NY10.inf <- subset(data, (data$strain == "NY10") & (data$cycle == FALSE)
                             & (data$infected == TRUE))

# Group infected and disseminated.
data.pO.WN02.inf <- aggregate(data.pO.WN02.inf[, c("ovi", "total")], by=list(data.pO.WN02.inf$temperature), FUN=sum)
colnames(data.pO.WN02.inf) <- c('temperature', 'ovi', 'total')
data.pO.WN02.inf$p.ovi <- data.pO.WN02.inf$ovi / data.pO.WN02.inf$total
data.pO.WN02.inf$se <- sqrt(data.pO.WN02.inf$p.ovi * (1 - data.pO.WN02.inf$p.ovi) / data.pO.WN02.inf$total)
data.pO.WN02.inf

data.pO.NY10.inf <- aggregate(data.pO.NY10.inf[, c("ovi", "total")], by=list(data.pO.NY10.inf$temperature), FUN=sum)
colnames(data.pO.NY10.inf) <- c('temperature', 'ovi', 'total')
data.pO.NY10.inf$p.ovi <- data.pO.NY10.inf$ovi / data.pO.NY10.inf$total
data.pO.NY10.inf$se <- sqrt(data.pO.NY10.inf$p.ovi * (1 - data.pO.NY10.inf$p.ovi) / data.pO.NY10.inf$total)
data.pO.NY10.inf

plot()

par(mfrow=c(2, 3))
errbar.length = 0.05
errbar.width = 1


plot(data.pO.unexp$temperature, data.pO.unexp$p.ovi, pch=20, ylim=c(0, 1),
      xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="Unexposed")
arrows(data.pO.unexp$temperature, data.pO.unexp$p.ovi, y1=data.pO.unexp$p.ovi + data.pO.unexp$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.unexp$temperature, data.pO.unexp$p.ovi, y1=data.pO.unexp$p.ovi - data.pO.unexp$se, 
       angle=90, length=errbar.length, lwd=errbar.width)


plot(data.pO.WN02.uninf$temperature, data.pO.WN02.uninf$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="WN02 (exp. not infected)")
arrows(data.pO.WN02.uninf$temperature, data.pO.WN02.uninf$p.ovi, y1=data.pO.WN02.uninf$p.ovi + data.pO.WN02.uninf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.WN02.uninf$temperature, data.pO.WN02.uninf$p.ovi, y1=data.pO.WN02.uninf$p.ovi - data.pO.WN02.uninf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)


plot(data.pO.NY10.uninf$temperature, data.pO.NY10.uninf$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="NY10 (exp. not infected)")
arrows(data.pO.NY10.uninf$temperature, data.pO.NY10.uninf$p.ovi, y1=data.pO.NY10.uninf$p.ovi + data.pO.NY10.uninf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.NY10.uninf$temperature, data.pO.NY10.uninf$p.ovi, y1=data.pO.NY10.uninf$p.ovi - data.pO.NY10.uninf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)


plot("")

plot(data.pO.WN02.inf$temperature, data.pO.WN02.inf$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="WN02 (infected)")
arrows(data.pO.WN02.inf$temperature, data.pO.WN02.inf$p.ovi, y1=data.pO.WN02.inf$p.ovi + data.pO.WN02.inf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.WN02.inf$temperature, data.pO.WN02.inf$p.ovi, y1=data.pO.WN02.inf$p.ovi - data.pO.WN02.inf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)

plot(data.pO.NY10.inf$temperature, data.pO.NY10.inf$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="NY10 (infected)")
arrows(data.pO.NY10.inf$temperature, data.pO.NY10.inf$p.ovi, y1=data.pO.NY10.inf$p.ovi + data.pO.NY10.inf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.NY10.inf$temperature, data.pO.NY10.inf$p.ovi, y1=data.pO.NY10.inf$p.ovi - data.pO.NY10.inf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)

sink("oviposition.txt")
cat("
    model{
    
    ## Priors
    Tmin ~ dnorm(8.2, 1/3^2)
    Tmax ~ dnorm(33.2, 1/3^2)
    rmax ~ dunif(0, 1)
    alpha ~ dunif(0, 1)
    beta ~ dgamma(0.35^2 / 0.2^2, 0.35 / 0.2^2)
    
    ## Derived Quantities and Predictions
    Topt <- alpha * Tmax + (1 - alpha) * Tmin
    s <- alpha * (1 - alpha) / beta^2
    
    ## Likelihood
    for(i in 1:N.obs) {
      p[i] <- max( (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (
                          - log(Tmax - Tmin) 
                        + alpha * log( max((temp[i] - Tmin) / alpha, 10^-20)) 
                        + (1 - alpha) * log( max((Tmax - temp[i]) / (1 - alpha), 10^-20) ))),
                        10^-10)
      n[i] ~ dbin(p[i], N[i])
    }
    } # close model
    ",fill=TRUE)
sink()

##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin = runif(1, min=0, max=10),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0, max=1),  
  beta = runif(1, min=0.1, max=0.9))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "Topt", "s")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains

##### Organize Data for JAGS        
jag.data <- list(n=data.pO.unexp$ovi, temp=data.pO.unexp$temperature,
                 N=data.pO.unexp$total,
                 N.obs=length(data.pO.unexp$temperature))

jag.data
##### Run JAGS
#unexp.pO.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                        model.file="oviposition.txt", n.thin=nt, n.chains=nc, 
#                        n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
unexp.pO.out <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                        model.file="oviposition.txt", 
                                        n.thin=nt, n.chains=nc, 
                                        n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                        jags.seed = 42))

unexp.pO.out
mcmcplot(unexp.pO.out)

##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin = runif(1, min=0, max=10),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0, max=1),  
  beta = runif(1, min=0.1, max=0.9))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "Topt", "s")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains

##### Organize Data for JAGS        
jag.data <- list(n=data.pO.WN02.uninf$ovi, temp=data.pO.WN02.uninf$temperature,
                 N=data.pO.WN02.uninf$total,
                 N.obs=length(data.pO.WN02.uninf$temperature))

jag.data

##### Run JAGS
#WN02.uninf.pO.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                     model.file="oviposition.txt", n.thin=nt, n.chains=nc, 
#                     n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
WN02.uninf.pO.out <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                            model.file="oviposition.txt", 
                                            n.thin=nt, n.chains=nc, 
                                            n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                            jags.seed = 42))
WN02.uninf.pO.out
mcmcplot(WN02.uninf.pO.out)


##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "Topt", "s")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains

##### Organize Data for JAGS        
jag.data <- list(n=data.pO.NY10.uninf$ovi, temp=data.pO.NY10.uninf$temperature,
                 N=data.pO.NY10.uninf$total,
                 N.obs=length(data.pO.NY10.uninf$temperature))

jag.data

##### Run JAGS
#NY10.uninf.pO.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                          model.file="oviposition.txt", n.thin=nt, n.chains=nc, 
#                          n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
NY10.uninf.pO.out <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                            model.file="oviposition.txt", 
                                            n.thin=nt, n.chains=nc, 
                                            n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                            jags.seed = 42))
NY10.uninf.pO.out
mcmcplot(NY10.uninf.pO.out)


##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin = runif(1, min=0, max=10),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0, max=1),  
  beta = runif(1, min=0.1, max=0.9))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "Topt", "s")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains

##### Organize Data for JAGS        
jag.data <- list(n=data.pO.WN02.inf$ovi, temp=data.pO.WN02.inf$temperature,
                 N=data.pO.WN02.inf$total,
                 N.obs=length(data.pO.WN02.inf$temperature))

jag.data

##### Run JAGS
#WN02.inf.pO.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                          model.file="oviposition.txt", n.thin=nt, n.chains=nc, 
#                          n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
WN02.inf.pO.out <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                                 model.file="oviposition.txt", 
                                                 n.thin=nt, n.chains=nc, 
                                                 n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                                 jags.seed = 42))

WN02.inf.pO.out
mcmcplot(WN02.inf.pO.out)


##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin = runif(1, min=0, max=10),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0, max=1),  
  beta = runif(1, min=0.1, max=0.9))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "Topt", "s")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains

##### Organize Data for JAGS        
jag.data <- list(n=data.pO.NY10.inf$ovi, temp=data.pO.NY10.inf$temperature,
                 N=data.pO.NY10.inf$total,
                 N.obs=length(data.pO.NY10.inf$temperature))

jag.data

##### Run JAGS
#NY10.inf.pO.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                        model.file="oviposition.txt", n.thin=nt, n.chains=nc, 
#                        n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
NY10.inf.pO.out <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                               model.file="oviposition.txt", 
                                               n.thin=nt, n.chains=nc, 
                                               n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                               jags.seed = 42))
NY10.inf.pO.out
mcmcplot(NY10.inf.pO.out)

par(mfrow=c(2, 3))

temps <- seq(0, 40, 0.1) 


plot(data.pO.unexp$temperature, data.pO.unexp$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="Unexposed")
arrows(data.pO.unexp$temperature, data.pO.unexp$p.ovi, y1=data.pO.unexp$p.ovi + data.pO.unexp$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.unexp$temperature, data.pO.unexp$p.ovi, y1=data.pO.unexp$p.ovi - data.pO.unexp$se, 
       angle=90, length=errbar.length, lwd=errbar.width)

chains.pO.unexp <- MCMCchains(unexp.pO.out, params=c("Tmin", "Tmax", "rmax", 
                                                     "alpha", "beta"))
curves <- apply(chains.pO.unexp, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))
# Find mean curve and credible intervals.
meancurve.pO.unexp <- apply(curves, 1, mean)
CI.pO.unexp <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.pO.unexp, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.pO.unexp [1, ],
                                rev(CI.pO.unexp [2, ])), 
        col=alpha("steelblue", 0.3), lty=0)


plot(data.pO.WN02.uninf$temperature, data.pO.WN02.uninf$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="WN02 (exp. not infected)")
arrows(data.pO.WN02.uninf$temperature, data.pO.WN02.uninf$p.ovi, y1=data.pO.WN02.uninf$p.ovi + data.pO.WN02.uninf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.WN02.uninf$temperature, data.pO.WN02.uninf$p.ovi, y1=data.pO.WN02.uninf$p.ovi - data.pO.WN02.uninf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)

chains.pO.WN02.uninf <- MCMCchains(WN02.uninf.pO.out, params=c("Tmin", "Tmax", "rmax", 
                                                     "alpha", "beta"))
curves <- apply(chains.pO.WN02.uninf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))
# Find mean curve and credible intervals.
meancurve.pO.WN02.uninf <- apply(curves, 1, mean)
CI.pO.WN02.uninf <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.pO.WN02.uninf, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.pO.WN02.uninf[1, ],
                                rev(CI.pO.WN02.uninf[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)

plot(data.pO.NY10.uninf$temperature, data.pO.NY10.uninf$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="NY10 (exp. not infected)")
arrows(data.pO.NY10.uninf$temperature, data.pO.NY10.uninf$p.ovi, y1=data.pO.NY10.uninf$p.ovi + data.pO.NY10.uninf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.NY10.uninf$temperature, data.pO.NY10.uninf$p.ovi, y1=data.pO.NY10.uninf$p.ovi - data.pO.NY10.uninf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)

chains.pO.NY10.uninf <- MCMCchains(NY10.uninf.pO.out, params=c("Tmin", "Tmax", "rmax", 
                                                               "alpha", "beta"))
curves <- apply(chains.pO.NY10.uninf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))
# Find mean curve and credible intervals.
meancurve.pO.NY10.uninf <- apply(curves, 1, mean)
CI.pO.NY10.uninf <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.pO.NY10.uninf, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.pO.NY10.uninf[1, ],
                                rev(CI.pO.NY10.uninf[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)

plot("")

plot(data.pO.WN02.inf$temperature, data.pO.WN02.inf$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="WN02 (infected)")
arrows(data.pO.WN02.inf$temperature, data.pO.WN02.inf$p.ovi, y1=data.pO.WN02.inf$p.ovi + data.pO.WN02.inf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.WN02.inf$temperature, data.pO.WN02.inf$p.ovi, y1=data.pO.WN02.inf$p.ovi - data.pO.WN02.inf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)

chains.pO.WN02.inf <- MCMCchains(WN02.inf.pO.out, params=c("Tmin", "Tmax", "rmax", 
                                                               "alpha", "beta"))
curves <- apply(chains.pO.WN02.inf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))
# Find mean curve and credible intervals.
meancurve.pO.WN02.inf <- apply(curves, 1, mean)
CI.pO.WN02.inf <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.pO.WN02.inf, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.pO.WN02.inf[1, ],
                                rev(CI.pO.WN02.inf[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)


plot(data.pO.NY10.inf$temperature, data.pO.NY10.inf$p.ovi, pch=20, ylim=c(0, 1),
     xlim=c(0, 40),
     xlab="Temperature [°C]", ylab="Proportion ovipositing",
     main="NY10 (infected)")
arrows(data.pO.NY10.inf$temperature, data.pO.NY10.inf$p.ovi, y1=data.pO.NY10.inf$p.ovi + data.pO.NY10.inf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.pO.NY10.inf$temperature, data.pO.NY10.inf$p.ovi, y1=data.pO.NY10.inf$p.ovi - data.pO.NY10.inf$se, 
       angle=90, length=errbar.length, lwd=errbar.width)

chains.pO.NY10.inf <- MCMCchains(NY10.inf.pO.out, params=c("Tmin", "Tmax", "rmax", 
                                                               "alpha", "beta"))
curves <- apply(chains.pO.NY10.inf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))
# Find mean curve and credible intervals.
meancurve.pO.NY10.inf <- apply(curves, 1, mean)
CI.pO.NY10.inf <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.pO.NY10.inf, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.pO.NY10.inf[1, ],
                                rev(CI.pO.NY10.inf[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)


par(mfrow=c(1,1))

plot(temps, meancurve.pO.unexp, col="steelblue", type='l', xlim=c(0, 40), 
     ylim=c(0, 1),
     xlab="Temperature [°C]", ylab="prop. ovipositing", lwd=2,
     main="Oviposition")
polygon(c(temps, rev(temps)), c(CI.pO.unexp [1, ],
                                rev(CI.pO.unexp [2, ])), 
        col=alpha("steelblue", 0.3), lty=0)

lines(temps, meancurve.pO.WN02.inf, col="darkgreen", lwd=2)
polygon(c(temps, rev(temps)), c(CI.pO.WN02.inf[1, ],
                                rev(CI.pO.WN02.inf[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)

lines(temps, meancurve.pO.NY10.inf, col="darkblue", lwd=2)
polygon(c(temps, rev(temps)), c(CI.pO.NY10.inf[1, ],
                                rev(CI.pO.NY10.inf[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)


saveRDS(unexp.pO.out, "./pO_unexposed.RDS")
saveRDS(WN02.uninf.pO.out, "./pO_WN02_uninfected.RDS")
saveRDS(NY10.uninf.pO.out, "./pO_NY10_uninfected.RDS")
saveRDS(WN02.inf.pO.out, "./pO_WN02_infected.RDS")
saveRDS(NY10.inf.pO.out, "./pO_NY10_infected.RDS")

unexp.pO.out <- readRDS("./pO_unexposed.RDS")
WN02.uninf.pO.out <- readRDS("./pO_WN02_uninfected.RDS")
NY10.uninf.pO.out <- readRDS("./pO_NY10_uninfected.RDS")
WN02.inf.pO.out <- readRDS("./pO_WN02_infected.RDS")
NY10.inf.pO.out <- readRDS("./pO_NY10_infected.RDS")

#View(unexp.pO.out$BUGSoutput$summary[c(2,3,1,7,4,5), c(1,2,3,7,8,9)])
#View(WN02.uninf.pO.out$BUGSoutput$summary[c(2,3,1,7,4,5), c(1,2,3,7,8,9)])
#View(WN02.inf.pO.out$BUGSoutput$summary[c(2,3,1,7,4,5), c(1,2,3,7,8,9)])
#View(NY10.uninf.pO.out$BUGSoutput$summary[c(2,3,1,7,4,5), c(1,2,3,7,8,9)])
#View(NY10.inf.pO.out$BUGSoutput$summary[c(2,3,1,7,4,5), c(1,2,3,7,8,9)])

get_param_estimates <- function(x) {
  return(x[c(2,3,1,7,4,5), c(1,3,7)])
}

par(mfrow=c(1, 2))
plot("", xlim=c(0, 40), ylim=c(0, 6), xlab="Temperature [°C]", ylab="",
     main="Cardinal and optimal temperatures",
      yaxt = "n")
axis(2, at=1:5, labels=c("unexp", "WN02n", "WN02i", "NY10n", "NY10i"), 
     las=1)
smry = list(unexp.pO.out$BUGSoutput$summary, WN02.uninf.pO.out$BUGSoutput$summary,
            WN02.inf.pO.out$BUGSoutput$summary,
            NY10.uninf.pO.out$BUGSoutput$summary,
            NY10.inf.pO.out$BUGSoutput$summary)
colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[1:3, 1], rep(i, 3), pch=20, col=colors[i])
  for(Tidx in 1:3) {
    lines(c(params[Tidx, 2], params[Tidx, 3]), rep(i, 2), lwd=2,
          col=alpha(colors[i], 0.5))
  }
}
eps = 0.2
text(10, 6 - eps, "Tmin")
text(24, 6 - eps, "Topt")
text(34, 6 - eps, "Tmax")

plot("", xlim=c(0, 1), ylim=c(0, 6), ylab="", xlab="prob. ovipositing",
     main="Peak oviposition",
     yaxt = "n")
axis(2, at=1:5, labels=c("unexp", "WN02n", "WN02i", "NY10n", "NY10i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[4, 1], i, pch=20, col=colors[i])
  lines(c(params[4, 2], params[4, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))

}

plot("", xlim=c(0, 6), ylim=c(0, 1), xlab="", ylab="alpha", xaxt = "n")
axis(1, at=1:5, labels=c("unexp", "WN02n", "WN02i", "NY10n", "NY10i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(i, params[5, 1], pch=20, col=colors[i])
  lines(rep(i, 2), c(params[5, 2], params[5, 3]), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

plot("", xlim=c(0, 6), ylim=c(0, 1), xlab="", ylab="beta",
     main="Proportion ovipositing", xaxt = "n")
axis(1, at=1:5, labels=c("unexp", "WN02n", "WN02i", "NY10n", "NY10i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(i, params[6, 1], pch=20, col=colors[i])
  lines(rep(i, 2), c(params[6, 2], params[6, 3]), lwd=2,
        col=alpha(colors[i], 0.5))
  
}
