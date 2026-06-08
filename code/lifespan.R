# Lifespan thermal response analysis
# Fits flexTPC models to mosquito lifespan data using JAGS.
# Run from the repository root directory.

set.seed(42) # Set seed for reproducibility.

# Load packages
library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

# FlexTPC model.
flexTPC <- function(T, Tmin, Tmax, rmax, alpha, beta) {
  s <- alpha * (1 - alpha) / beta^2
  result <- rep(0, length(T))
  Tidx = (T > Tmin) & (T < Tmax)
  result[Tidx] <- rmax * exp(s * (alpha * log( (T[Tidx] - Tmin) / alpha) 
                                  + (1 - alpha) * log( (Tmax - T[Tidx]) / (1 - alpha))
                                  - log(Tmax - Tmin)) ) 
  return(result)
}

data.lf <- read.csv("./data/lifespan.csv")
data.lf

summary(data.lf$lifespan)

# Unexposed lifespan, constant temperature
data.lf.unexp <- subset(data.lf, ((data.lf$genotype == "none") &
                                    (data.lf$cycle == FALSE) &
                                    (data.lf$infected == FALSE)))

# Uninfected lifespan, constant temperature
data.lf.WN02.uninf <- subset(data.lf, ((data.lf$genotype == "WN02") &
                                         (data.lf$cycle == FALSE) &
                                         (data.lf$infected == FALSE)))
data.lf.NY10.uninf <- subset(data.lf, ((data.lf$genotype == "NY10") &
                                         (data.lf$cycle == FALSE) &
                                         (data.lf$infected == FALSE)))

# Infected lifespan, constant temperature
data.lf.WN02.inf <- subset(data.lf, ((data.lf$genotype == "WN02") &
                                       (data.lf$cycle == FALSE) &
                                       (data.lf$infected == TRUE)))
data.lf.NY10.inf <- subset(data.lf, ((data.lf$genotype == "NY10") &
                                       (data.lf$cycle == FALSE) &
                                       (data.lf$infected == TRUE)))

par(mfrow=c(2, 3))
plot(data.lf.unexp$temperature, data.lf.unexp$lifespan,
     col=alpha("black", 0.3), xlim=c(0, 40), ylim=c(0, 100),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="Unexposed")
plot(data.lf.WN02.uninf$temperature, data.lf.WN02.uninf$lifespan,
     col=alpha("black", 0.3), xlim=c(0, 40), ylim=c(0, 100),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="WN02 (never infected)")
plot(data.lf.NY10.uninf$temperature, data.lf.NY10.uninf$lifespan,
     col=alpha("black", 0.3), xlim=c(0, 40), ylim=c(0, 100),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="NY10 (never infected)")

plot("")
plot(data.lf.WN02.inf$temperature, data.lf.WN02.inf$lifespan,
     col=alpha("black", 0.3), xlim=c(0, 40), ylim=c(0, 100),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="WN02 (infected+disseminated)")
plot(data.lf.NY10.inf$temperature, data.lf.NY10.inf$lifespan,
     col=alpha("black", 0.3), xlim=c(0, 40), ylim=c(0, 100),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="NY10 (infected+disseminated)")


par(mfrow=c(2, 3))
errbar.length=0.06
errbar.width=2
temperatures <- unique(data.lf$temperature)
temperatures

m <- with(data.lf.unexp, tapply(lifespan, temperature, mean))
s <- with(data.lf.unexp, tapply(lifespan, temperature, sd))
n <- with(data.lf.unexp, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="Unexposed", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

m <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, mean))
s <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, sd))
n <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="WN02 (never infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

m <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, mean))
s <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, sd))
n <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="NY10 (never infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

plot("")

m <- with(data.lf.WN02.inf, tapply(lifespan, temperature, mean))
s <- with(data.lf.WN02.inf, tapply(lifespan, temperature, sd))
n <- with(data.lf.WN02.inf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="WN02 (infected+disseminated)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

m <- with(data.lf.NY10.inf, tapply(lifespan, temperature, mean))
s <- with(data.lf.NY10.inf, tapply(lifespan, temperature, sd))
n <- with(data.lf.NY10.inf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="NY10 (infected+disseminated)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)


sink("lifespan_nb.txt")
cat("
    model{
    
    ## Priors
    rmax ~ dunif(0, 150)
    
    Tmin ~ dnorm(5, 1/2.5^2)
    Tmax ~ dnorm(34.9, 1 /2.5^2) # From Shocket et al.
    
    alpha ~ dunif(0, 1)
    Topt <- alpha * Tmax + (1 - alpha) * Tmin
    
    beta ~ dgamma(0.35^2 / 0.2^2, 0.35 / 0.2^2)
    
    s <- alpha * (1 - alpha) / beta^2
    
    # Size parameter
    r ~ dunif(0, 50)
    
    ## Likelihood
    for(i in 1:N.obs){
    trait.mu[i] <-  max( (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (
                        - log(Tmax - Tmin) 
                      + alpha * log( max((temp[i] - Tmin) / alpha, 10^-15)) 
                      + (1 - alpha) * log( max((Tmax - temp[i]) / (1 - alpha), 10^-15) ))),
                      10^-6)
    p[i] <- r/(r + trait.mu[i])
    trait[i] ~ dnegbin(p[i], r)
    }
    
    } # close model
    ",fill=TRUE)
sink()

data.lf.unexp
for(temp in temperatures) {
  hist(data.lf.unexp$lifespan[data.lf.unexp$temperature == temp],
       main=paste("Unexposed ", temp, "°C", sep=""), xlab="Temperature [°C]",
       xlim=c(0, 80))
}

for(temp in temperatures) {
  hist(data.lf.WN02.uninf$lifespan[data.lf.WN02.uninf$temperature == temp],
       main=paste("WN02 NI ", temp, "°C", sep=""), xlab="Temperature [°C]",
       xlim=c(0, 80))
}

for(temp in temperatures) {
  hist(data.lf.NY10.uninf$lifespan[data.lf.NY10.uninf$temperature == temp],
       main=paste("NY10 NI ", temp, "°C", sep=""), xlab="Temperature [°C]",
       xlim=c(0, 80))
}


for(temp in temperatures) {
  hist(data.lf.WN02.inf$lifespan[data.lf.WN02.inf$temperature == temp],
       main=paste("WN02 INF ", temp, "°C", sep=""), xlab="Temperature [°C]",
       xlim=c(0, 80))
}

for(temp in temperatures) {
  hist(data.lf.NY10.inf$lifespan[data.lf.NY10.inf$temperature == temp],
       main=paste("NY10 INF ", temp, "°C", sep=""), xlab="Temperature [°C]",
       xlim=c(0, 80))
}

### Unexposed lifespan

##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin = runif(1, min=0, max=5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=50),
  alpha = runif(1, min=0.3, max=0.7),  
  beta = runif(1, min=0.2, max=0.5),
  r = runif(1, 1, 50))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "r", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains


##### Organize Data for JAGS        
jag.data <- list(trait=data.lf.unexp$lifespan, temp=data.lf.unexp$temperature,
                 N.obs=length(data.lf.unexp$temperature))

data.lf.unexp$lifespan
##### Run JAGS
unexp.lf.out.nb <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                               model.file="lifespan_nb.txt", n.thin=nt, n.chains=nc, 
                                               n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                               jags.seed = 42))
unexp.lf.out.nb
mcmcplot(unexp.lf.out.nb)


chains.unexp.lf <- MCMCchains(unexp.lf.out.nb, params=c("Tmin", "Tmax", "rmax", 
                                                        "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.unexp.lf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.unexp.lf <- apply(curves, 1, mean)
CI.unexp.lf <- apply(curves, 1, quantile, c(0.025, 0.975))

par(mfrow=c(1, 1))

m <- with(data.lf.unexp, tapply(lifespan, temperature, mean))
s <- with(data.lf.unexp, tapply(lifespan, temperature, sd))
n <- with(data.lf.unexp, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="Unexposed", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.unexp.lf, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.unexp.lf[1, ],
                                rev(CI.unexp.lf[2,])), 
        col=alpha("steelblue", 0.3), lty=0)


##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin = runif(1, min=0, max=5),
  Tmax = runif(1, min=35, max=40),
  rmax = runif(1, min=0, max=50),
  alpha = runif(1, min=0.3, max=0.7),  
  beta = runif(1, min=0.2, max=0.5),
  r = runif(1, 1, 50))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "r", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains


##### Organize Data for JAGS        
jag.data <- list(trait=data.lf.WN02.uninf$lifespan, temp=data.lf.WN02.uninf$temperature,
                 N.obs=length(data.lf.WN02.uninf$temperature))

data.lf.WN02.uninf$lifespan
##### Run JAGS
WN02.uninf.lf.nb <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                                model.file="lifespan_nb.txt", n.thin=nt, n.chains=nc, 
                                                n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                                jags.seed = 42))
WN02.uninf.lf.nb
mcmcplot(WN02.uninf.lf.nb)


chains.WN02.uninf.lf.nb <- MCMCchains(WN02.uninf.lf.nb, params=c("Tmin", "Tmax", "rmax", 
                                                                 "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.WN02.uninf.lf.nb, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.WN02.uninf.lf.nb <- apply(curves, 1, mean)
CI.WN02.uninf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))

par(mfrow=c(1, 1))

m <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, mean))
s <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, sd))
n <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="WN02 (never infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.uninf.lf.nb, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.WN02.uninf.lf.nb[1, ],
                                rev(CI.WN02.uninf.lf.nb[2, ])), 
        col=alpha("steelblue", 0.3), lty=0)


## ##### Organize Data for JAGS        
jag.data <- list(trait=data.lf.NY10.uninf$lifespan, temp=data.lf.NY10.uninf$temperature,
                 N.obs=length(data.lf.NY10.uninf$temperature))
data.lf.NY10.uninf

data.lf.unexp$lifespan
##### Run JAGS
NY10.uninf.lf.nb <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                                model.file="lifespan_nb.txt", n.thin=nt, n.chains=nc, 
                                                n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                                jags.seed = 42))
NY10.uninf.lf.nb
mcmcplot(NY10.uninf.lf.nb)


chains.NY10.uninf.lf.nb <- MCMCchains(NY10.uninf.lf.nb, params=c("Tmin", "Tmax", "rmax", 
                                                                 "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.NY10.uninf.lf.nb, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.NY10.uninf.lf.nb <- apply(curves, 1, mean)
CI.NY10.uninf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))

par(mfrow=c(1, 1))

m <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, mean))
s <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, sd))
n <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="NY10 (never infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.uninf.lf.nb, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.NY10.uninf.lf.nb[1, ],
                                rev(CI.NY10.uninf.lf.nb[2, ])), 
        col=alpha("steelblue", 0.3), lty=0)

## ##### Organize Data for JAGS        
jag.data <- list(trait=data.lf.WN02.inf$lifespan, temp=data.lf.WN02.inf$temperature,
                 N.obs=length(data.lf.WN02.inf$temperature))

##### Run JAGS
WN02.inf.lf.nb <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                              model.file="lifespan_nb.txt", n.thin=nt, n.chains=nc, 
                                              n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                              jags.seed = 42))
WN02.inf.lf.nb
mcmcplot(WN02.inf.lf.nb)


chains.WN02.inf.lf.nb <- MCMCchains(WN02.inf.lf.nb, params=c("Tmin", "Tmax", "rmax", 
                                                             "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.WN02.inf.lf.nb, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.WN02.inf.lf.nb <- apply(curves, 1, mean)
CI.WN02.inf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))

par(mfrow=c(1, 1))

m <- with(data.lf.WN02.inf, tapply(lifespan, temperature, mean))
s <- with(data.lf.WN02.inf, tapply(lifespan, temperature, sd))
n <- with(data.lf.WN02.inf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="WN02 (infected+disseminated)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.inf.lf.nb, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.WN02.inf.lf.nb[1, ],
                                rev(CI.WN02.inf.lf.nb[2, ])), 
        col=alpha("steelblue", 0.3), lty=0)

## ##### Organize Data for JAGS        
jag.data <- list(trait=data.lf.NY10.inf$lifespan, temp=data.lf.NY10.inf$temperature,
                 N.obs=length(data.lf.NY10.inf$temperature))

##### Run JAGS
NY10.inf.lf.nb <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                              model.file="lifespan_nb.txt", n.thin=nt, n.chains=nc, 
                                              n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                              jags.seed = 42))
NY10.inf.lf.nb
mcmcplot(NY10.inf.lf.nb)


chains.NY10.inf.lf.nb <- MCMCchains(NY10.inf.lf.nb, params=c("Tmin", "Tmax", "rmax", 
                                                             "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.NY10.inf.lf.nb, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.NY10.inf.lf.nb <- apply(curves, 1, mean)
CI.NY10.inf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))

par(mfrow=c(1, 1))

m <- with(data.lf.NY10.inf, tapply(lifespan, temperature, mean))
s <- with(data.lf.NY10.inf, tapply(lifespan, temperature, sd))
n <- with(data.lf.NY10.inf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="NY10 (infected+disseminated)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.inf.lf.nb, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.NY10.inf.lf.nb[1, ],
                                rev(CI.NY10.inf.lf.nb[2, ])), 
        col=alpha("steelblue", 0.3), lty=0)



par(mfrow=c(2, 3))
errbar.length=0.06
errbar.width=2
temperatures <- unique(data.lf$temperature)
temperatures

m <- with(data.lf.unexp, tapply(lifespan, temperature, mean))
s <- with(data.lf.unexp, tapply(lifespan, temperature, sd))
n <- with(data.lf.unexp, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="Unexposed", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)


# Plot mean curve and 95% credible interval.
lines(temps, meancurve.unexp.lf, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.unexp.lf[1, ],
                                rev(CI.unexp.lf[2,])), 
        col=alpha("steelblue", 0.3), lty=0)


m <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, mean))
s <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, sd))
n <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="WN02 (never infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.uninf.lf.nb, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.uninf.lf.nb[1, ],
                                rev(CI.WN02.uninf.lf.nb[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)

m <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, mean))
s <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, sd))
n <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="NY10 (never infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

lines(temps, meancurve.NY10.uninf.lf.nb, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.uninf.lf.nb[1, ],
                                rev(CI.NY10.uninf.lf.nb[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)

plot("")

m <- with(data.lf.WN02.inf, tapply(lifespan, temperature, mean))
s <- with(data.lf.WN02.inf, tapply(lifespan, temperature, sd))
n <- with(data.lf.WN02.inf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="WN02 (infected+disseminated)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

lines(temps, meancurve.WN02.inf.lf.nb, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.inf.lf.nb[1, ],
                                rev(CI.WN02.inf.lf.nb[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)


m <- with(data.lf.NY10.inf, tapply(lifespan, temperature, mean))
s <- with(data.lf.NY10.inf, tapply(lifespan, temperature, sd))
n <- with(data.lf.NY10.inf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="NY10 (infected+disseminated)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

lines(temps, meancurve.NY10.inf.lf.nb, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.inf.lf.nb[1, ],
                                rev(CI.NY10.inf.lf.nb[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)


par(mfrow=c(1,1))

plot(temps, meancurve.unexp.lf, col="steelblue", type='l', xlim=c(0, 40), 
     ylim=c(0, 30),
     xlab="Temperature [°C]", ylab="Lifespan [days]", lwd=2,
     main="Unexposed vs infected")
polygon(c(temps, rev(temps)), c(CI.unexp.lf[1, ],
                                rev(CI.unexp.lf[2, ])), 
        col=alpha("steelblue", 0.2), lty=0)

lines(temps, meancurve.WN02.inf.lf.nb, col="darkgreen", lwd=2)
polygon(c(temps, rev(temps)), c(CI.WN02.inf.lf.nb[1, ],
                                rev(CI.WN02.inf.lf.nb[2, ])), 
        col=alpha("darkgreen", 0.2), lty=0)

lines(temps, meancurve.NY10.inf.lf.nb, col="darkblue", lwd=2)
polygon(c(temps, rev(temps)), c(CI.NY10.inf.lf.nb[1, ],
                                rev(CI.NY10.inf.lf.nb[2, ])), 
        col=alpha("darkblue", 0.2), lty=0)


saveRDS(unexp.lf.out.nb, "./lf_unexposed.RDS")
saveRDS(WN02.uninf.lf.nb, "./lf_WN02_uninfected.RDS")
saveRDS(NY10.uninf.lf.nb, "./lf_NY10_uninfected.RDS")
saveRDS(WN02.inf.lf.nb, "./lf_WN02_infected.RDS")
saveRDS(NY10.inf.lf.nb, "./lf_NY10_infected.RDS")

NY10.inf.lf.nb