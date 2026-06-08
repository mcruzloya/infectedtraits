# Bloodfeeding thermal response analysis
# Summarizes bloodfeeding data by temperature and infection status,
# fits flexTPC models using JAGS, and saves posterior model fits.
# Run from the repository root directory.

set.seed(42) # Set seed for reproducibility.

# Import packages.
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

data.bf <- read.csv("./data/bloodfeeding_lf.csv")
data.bf$infected <- as.logical(data.bf$infected)
data.bf$disseminated <- as.logical(data.bf$disseminated)
data.bf$cycle <- as.logical(data.bf$cycle)

# Unexposed BF, constant temperature
data.bf.unexp <- subset(data.bf, ((data.bf$strain == "None") &
                                    (data.bf$cycle == FALSE) &
                                    (data.bf$infected == FALSE)))


temperatures <- c(10, 15, 20, 25, 30, 33)
data.bf.unexp.smry <- data.frame(temp=temperatures, bloodmeals=NaN, 
                                bloodmeals_offered=NaN)
data.bf.unexp.smry
for(temp in temperatures) {
  df <- subset(data.bf.unexp, data.bf.unexp$temperature == temp)
  data.bf.unexp.smry$bloodmeals[data.bf.unexp.smry$temp == temp]  <- sum(df$bloodmeals)
  data.bf.unexp.smry$bloodmeals_offered[data.bf.unexp.smry$temp == temp]  <- sum(df$bloodmeals_offered)
}
data.bf.unexp.smry$prop <- data.bf.unexp.smry$bloodmeals / data.bf.unexp.smry$bloodmeals_offered
data.bf.unexp.smry$sem <- sqrt(data.bf.unexp.smry$prop * (1 - data.bf.unexp.smry$prop) / data.bf.unexp.smry$bloodmeals_offered)


errbar.length=0.06
errbar.width=2

plot(temperatures, data.bf.unexp.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab='bites / (mosquito * week)', ylim=c(0, 0.4))
arrows(temperatures, data.bf.unexp.smry$prop, y1= data.bf.unexp.smry$prop + data.bf.unexp.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.unexp.smry$prop, y1= data.bf.unexp.smry$prop - data.bf.unexp.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Uninfected BF, constant temperature
data.bf.WN02.uninf <- subset(data.bf, ((data.bf$strain == "WN02") &
                                   (data.bf$cycle == FALSE) &
                                   (data.bf$infected == FALSE)))
data.bf.WN02.uninf.smry <- data.frame(temp=temperatures, bloodmeals=NaN, 
                                      bloodmeals_offered=NaN)
for(temp in temperatures) {
  df <- subset(data.bf.WN02.uninf, data.bf.WN02.uninf$temperature == temp)
  data.bf.WN02.uninf.smry$bloodmeals[data.bf.WN02.uninf.smry$temp == temp]  <- sum(df$bloodmeals)
  data.bf.WN02.uninf.smry$bloodmeals_offered[data.bf.WN02.uninf.smry$temp == temp]  <- sum(df$bloodmeals_offered)
}
data.bf.WN02.uninf.smry$prop <- data.bf.WN02.uninf.smry$bloodmeals / data.bf.WN02.uninf.smry$bloodmeals_offered
data.bf.WN02.uninf.smry$sem <- sqrt(data.bf.WN02.uninf.smry$prop * (1 - data.bf.WN02.uninf.smry$prop) / data.bf.WN02.uninf.smry$bloodmeals_offered)

data.bf.NY10.uninf <- subset(data.bf, ((data.bf$strain == "NY10") &
                                         (data.bf$cycle == FALSE) &
                                         (data.bf$infected == FALSE)))
data.bf.NY10.uninf.smry <- data.frame(temp=temperatures, bloodmeals=NaN, 
                                      bloodmeals_offered=NaN)
for(temp in temperatures) {
  df <- subset(data.bf.NY10.uninf, data.bf.NY10.uninf$temperature == temp)
  data.bf.NY10.uninf.smry$bloodmeals[data.bf.NY10.uninf.smry$temp == temp]  <- sum(df$bloodmeals)
  data.bf.NY10.uninf.smry$bloodmeals_offered[data.bf.NY10.uninf.smry$temp == temp]  <- sum(df$bloodmeals_offered)
}
data.bf.NY10.uninf.smry$prop <- data.bf.NY10.uninf.smry$bloodmeals / data.bf.NY10.uninf.smry$bloodmeals_offered
data.bf.NY10.uninf.smry$sem <- sqrt(data.bf.NY10.uninf.smry$prop * (1 - data.bf.NY10.uninf.smry$prop) / data.bf.NY10.uninf.smry$bloodmeals_offered)


# Infected bloodfeeding, constant temperature
data.bf.WN02.inf <- subset(data.bf, ((data.bf$strain == "WN02") &
                                         (data.bf$cycle == FALSE) &
                                         (data.bf$infected == TRUE)))

data.bf.WN02.inf.smry <- data.frame(temp=temperatures, bloodmeals=NaN, 
                                      bloodmeals_offered=NaN)
for(temp in temperatures) {
  df <- subset(data.bf.WN02.inf, data.bf.WN02.inf$temperature == temp)
  data.bf.WN02.inf.smry$bloodmeals[data.bf.WN02.inf.smry$temp == temp]  <- sum(df$bloodmeals)
  data.bf.WN02.inf.smry$bloodmeals_offered[data.bf.WN02.inf.smry$temp == temp]  <- sum(df$bloodmeals_offered)
}
data.bf.WN02.inf.smry$prop <- data.bf.WN02.inf.smry$bloodmeals / data.bf.WN02.inf.smry$bloodmeals_offered
data.bf.WN02.inf.smry$sem <- sqrt(data.bf.WN02.inf.smry$prop * (1 - data.bf.WN02.inf.smry$prop) / data.bf.WN02.inf.smry$bloodmeals_offered)
data.bf.WN02.inf.smry


data.bf.NY10.inf <- subset(data.bf, ((data.bf$strain == "NY10") &
                                         (data.bf$cycle == FALSE) &
                                         (data.bf$infected == TRUE)))
data.bf.NY10.inf.smry <- data.frame(temp=temperatures, bloodmeals=NaN, 
                                    bloodmeals_offered=NaN)
for(temp in temperatures) {
  df <- subset(data.bf.NY10.inf, data.bf.NY10.inf$temperature == temp)
  data.bf.NY10.inf.smry$bloodmeals[data.bf.NY10.inf.smry$temp == temp]  <- sum(df$bloodmeals)
  data.bf.NY10.inf.smry$bloodmeals_offered[data.bf.NY10.inf.smry$temp == temp]  <- sum(df$bloodmeals_offered)
}
data.bf.NY10.inf.smry$prop <- data.bf.NY10.inf.smry$bloodmeals / data.bf.NY10.inf.smry$bloodmeals_offered
data.bf.NY10.inf.smry$sem <- sqrt(data.bf.NY10.inf.smry$prop * (1 - data.bf.NY10.inf.smry$prop) / data.bf.NY10.inf.smry$bloodmeals_offered)
data.bf.NY10.inf.smry



par(mfrow=c(2, 3))
plot(temperatures, data.bf.unexp.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4),
     main="Unexposed")
arrows(temperatures, data.bf.unexp.smry$prop, y1= data.bf.unexp.smry$prop + data.bf.unexp.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.unexp.smry$prop, y1= data.bf.unexp.smry$prop - data.bf.unexp.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

plot(temperatures, data.bf.WN02.uninf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4),
     main="WN02 (uninfected)")
arrows(temperatures, data.bf.WN02.uninf.smry$prop, y1= data.bf.WN02.uninf.smry$prop + data.bf.WN02.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.WN02.uninf.smry$prop, y1= data.bf.WN02.uninf.smry$prop - data.bf.WN02.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

plot(temperatures, data.bf.NY10.uninf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4),
     main="NY10 (uninfected)")
arrows(temperatures, data.bf.NY10.uninf.smry$prop, y1= data.bf.NY10.uninf.smry$prop + data.bf.NY10.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.NY10.uninf.smry$prop, y1= data.bf.NY10.uninf.smry$prop - data.bf.NY10.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)



plot("")
plot(temperatures, data.bf.WN02.inf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4),
     main="WN02 (infected)")
arrows(temperatures, data.bf.WN02.inf.smry$prop, y1= data.bf.WN02.inf.smry$prop + data.bf.WN02.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.WN02.inf.smry$prop, y1= data.bf.WN02.inf.smry$prop - data.bf.WN02.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

plot(temperatures, data.bf.NY10.inf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4),
     main="NY10 (infected)")
arrows(temperatures, data.bf.NY10.inf.smry$prop, y1= data.bf.NY10.inf.smry$prop + data.bf.NY10.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.NY10.inf.smry$prop, y1= data.bf.NY10.inf.smry$prop - data.bf.NY10.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)


data.bf.unexp.smry
data.bf.WN02.uninf.smry
data.bf.WN02.inf.smry
data.bf.NY10.uninf.smry
data.bf.NY10.inf.smry

#write.csv(data.bf.unexp.smry, file="bf_prop_unexp.csv")
#write.csv(data.bf.WN02.uninf.smry, file="bf_prop_WN02_uninf.csv")
#write.csv(data.bf.WN02.inf.smry, file="bf_prop_WN02_inf.csv")
#write.csv(data.bf.NY10.uninf.smry, file="bf_prop_NY10_uninf.csv")
#write.csv(data.bf.NY10.inf.smry, file="bf_prop_NY10_inf.csv")

sink("bloodfeeding.txt")
cat("
    model{
    
    ## Priors
    rmax ~ dunif(0, 1)
    
    Tmin ~ dnorm(7.5, 1 / 3.75^2)
    Tmax ~ dnorm(35, 1 / 3^2) 
    alpha ~ dunif(0, 1)
    Topt <- alpha * Tmax + (1 - alpha) * Tmin
    
    beta ~ dgamma(0.35^2 / 0.2^2, 0.35 / 0.2^2)
    
    s <- alpha * (1 - alpha) / beta^2

    ## Likelihood
    for(i in 1:N.obs){
      p[i] <-  max( (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (
                          - log(Tmax - Tmin) 
                        + alpha * log( max((temp[i] - Tmin) / alpha, 10^-20)) 
                        + (1 - alpha) * log( max((Tmax - temp[i]) / (1 - alpha), 10^-20) ))),
                        10^-20)
      n[i] ~ dbinom(p[i], N[i])
    }
    
    } # close model
    ",fill=TRUE)
sink()


sink("bloodfeeding_unexp_priors.txt")
cat("
    model{
    
    ## Priors
    rmax ~ dunif(0, 1)
    
    Tmin ~ dnorm(9.385, 1 / 3.5^2)
    Tmax ~ dnorm(34.732, 1 / 3^2) 
    
    kappa <- 0.512 * (1 - 0.512) / 0.15^2 - 1.0
    alpha ~ dbeta(0.512 * kappa, (1.0 - 0.512) * kappa)
    Topt <- alpha * Tmax + (1 - alpha) * Tmin
    
    beta ~ dgamma(0.209^2 / 0.062^2, 0.209 / 0.062^2)
    
    s <- alpha * (1 - alpha) / beta^2

    ## Likelihood
    for(i in 1:N.obs){
      p[i] <-  max( (Tmax > temp[i]) * (Tmin < temp[i]) * rmax * exp(s * (
                          - log(Tmax - Tmin) 
                        + alpha * log( max((temp[i] - Tmin) / alpha, 10^-20)) 
                        + (1 - alpha) * log( max((Tmax - temp[i]) / (1 - alpha), 10^-20) ))),
                        10^-20)
      n[i] ~ dbinom(p[i], N[i])
    }
    
    } # close model
    ",fill=TRUE)
sink()


### Unexposed

##### Calculate initial values for MCMC.
## We are picking random values so that every chain will start at a different place. 
inits<-function(){list(
  Tmin = runif(1, min=5, max=10),
  Tmax = runif(1, min=33, max=35),
  rmax = runif(1, min=0, max=1),
  alpha = runif(1, min=0, max=1),  
  beta = runif(1, min=0.1, max=0.5))}

##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "Topt", "s")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains

data.bf.unexp.smry$bloodmeals
data.bf.unexp.smry$bloodmeals_offered
data.bf.unexp.smry$temp
##### Organize Data for JAGS        
jag.data <- list(temp=data.bf.unexp.smry$temp, n=data.bf.unexp.smry$bloodmeals, 
                 N=data.bf.unexp.smry$bloodmeals_offered,
                 N.obs=length(data.bf.unexp.smry$temp))

##### Run JAGS
#unexp.bf.out <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                     model.file="bloodfeeding.txt", n.thin=nt, n.chains=nc, 
#                     n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())

unexp.bf.out <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                               model.file="bloodfeeding.txt", n.thin=nt, n.chains=nc, 
                                               n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                               jags.seed = 42))
unexp.bf.out
#View(unexp.bf.out$BUGSoutput$summary[c("Tmin", "Topt", "Tmax", "rmax", "alpha", "beta"), 
#                  ])
mcmcplot(unexp.bf.out)


chains.unexp.bf <- MCMCchains(unexp.bf.out, params=c("Tmin", "Tmax", "rmax", 
                                                     "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.unexp.bf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                            x[5]))

# Find mean curve and credible intervals.
meancurve.unexp.bf <- apply(curves, 1, mean)
CI.unexp.bf <- apply(curves, 1, quantile, c(0.025, 0.975))

#par(mfrow=c(1, 1))

plot(temperatures, data.bf.unexp.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", xlim=c(0, 40), ylim=c(0, 0.4),
     main="Unexposed")
arrows(temperatures, data.bf.unexp.smry$prop, y1= data.bf.unexp.smry$prop + data.bf.unexp.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.unexp.smry$prop, y1= data.bf.unexp.smry$prop - data.bf.unexp.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.unexp.bf, col="steelblue")
polygon(c(temps, rev(temps)), c(CI.unexp.bf[1, ],
                                rev(CI.unexp.bf[2,])), 
        col=alpha("steelblue", 0.3), lty=0)

### WN02 unexposed


##### Parameters to Estimate
parameters <- c("Tmin", "Tmax", "rmax", "alpha", "beta", "s", "Topt")

##### MCMC Settings
# Number of posterior dist elements = [(ni - nb) / nt ] * nc = [ (25000 - 5000) / 8 ] * 3 = 7500
ni <- 300000 # number of iterations in each chain
nb <- 50000 # number of 'burn in' iterations to discard
nt <- 8 # thinning rate - jags saves every nt iterations in each chain
nc <- 8 # number of chains


##### Organize Data for JAGS        
jag.data <- list(n=data.bf.WN02.uninf.smry$bloodmeals, temp=data.bf.WN02.uninf.smry$temp,
                 N=data.bf.WN02.uninf.smry$bloodmeals_offered,
                 N.obs=length(data.bf.WN02.uninf.smry$temp))

jag.data

##### Run JAGS
#WN02.uninf.bf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                     model.file="bloodfeeding_unexp_priors.txt", n.thin=nt, n.chains=nc, 
#                     n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
WN02.uninf.bf <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                            model.file="bloodfeeding_unexp_priors.txt", 
                                            n.thin=nt, n.chains=nc, 
                                            n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                            jags.seed = 42))

WN02.uninf.bf
mcmcplot(WN02.uninf.bf)


chains.WN02.uninf.bf <- MCMCchains(WN02.uninf.bf, params=c("Tmin", "Tmax", "rmax", 
                                                     "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.WN02.uninf.bf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.WN02.uninf.bf <- apply(curves, 1, mean)
CI.WN02.uninf.bf <- apply(curves, 1, quantile, c(0.025, 0.975))




plot(temperatures, data.bf.WN02.uninf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4), xlim=c(0, 40),
     main="WN02 (uninfected)")
arrows(temperatures, data.bf.WN02.uninf.smry$prop, y1= data.bf.WN02.uninf.smry$prop + data.bf.WN02.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.WN02.uninf.smry$prop, y1= data.bf.WN02.uninf.smry$prop - data.bf.WN02.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)


jag.data

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.uninf.bf, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.uninf.bf[1, ],
                              rev(CI.WN02.uninf.bf[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)


## ##### Organize Data for JAGS        
jag.data <- list(n=data.bf.NY10.uninf.smry$bloodmeals, temp=data.bf.NY10.uninf.smry$temp,
                 N=data.bf.NY10.uninf.smry$bloodmeals_offered,
                 N.obs=length(data.bf.NY10.uninf.smry$temp))
jag.data

##### Run JAGS
#NY10.uninf.bf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                    model.file="bloodfeeding_unexp_priors.txt", n.thin=nt, n.chains=nc, 
#                    n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
NY10.uninf.bf <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                             model.file="bloodfeeding_unexp_priors.txt", 
                                             n.thin=nt, n.chains=nc, 
                                             n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                             jags.seed = 42))
NY10.uninf.bf
mcmcplot(NY10.uninf.bf)


chains.NY10.uninf.bf <- MCMCchains(NY10.uninf.bf, params=c("Tmin", "Tmax", "rmax", 
                                                       "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.NY10.uninf.bf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.NY10.uninf.bf <- apply(curves, 1, mean)
CI.NY10.uninf.bf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot(temperatures, data.bf.NY10.uninf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4), xlim=c(0, 40),
     main="NY10 (uninfected)")
arrows(temperatures, data.bf.NY10.uninf.smry$prop, y1= data.bf.NY10.uninf.smry$prop + data.bf.NY10.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.NY10.uninf.smry$prop, y1= data.bf.NY10.uninf.smry$prop - data.bf.NY10.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)


# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.uninf.bf, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.uninf.bf[1, ],
                                rev(CI.NY10.uninf.bf[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)

## ##### Organize Data for JAGS        
jag.data <- list(n=data.bf.WN02.inf.smry$bloodmeals, temp=data.bf.WN02.inf.smry$temp,
                 N=data.bf.WN02.inf.smry$bloodmeals_offered,
                 N.obs=length(data.bf.WN02.inf.smry$temp))
jag.data



##### Run JAGS
#WN02.inf.bf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                       model.file="bloodfeeding_unexp_priors.txt", n.thin=nt, n.chains=nc, 
#                       n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
WN02.inf.bf <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                             model.file="bloodfeeding_unexp_priors.txt", 
                                             n.thin=nt, n.chains=nc, 
                                             n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                             jags.seed = 42))

WN02.inf.bf

mcmcplot(WN02.inf.bf)


chains.WN02.inf.bf <- MCMCchains(WN02.inf.bf, params=c("Tmin", "Tmax", "rmax", 
                                                          "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.WN02.inf.bf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.WN02.inf.bf <- apply(curves, 1, mean)
CI.WN02.inf.bf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot("")
plot(temperatures, data.bf.WN02.inf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4), xlim=c(0, 40),
     main="WN02 (infected)")
arrows(temperatures, data.bf.WN02.inf.smry$prop, y1= data.bf.WN02.inf.smry$prop + data.bf.WN02.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.WN02.inf.smry$prop, y1= data.bf.WN02.inf.smry$prop - data.bf.WN02.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.inf.bf, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.inf.bf[1, ],
                                rev(CI.WN02.inf.bf[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)

## ##### Organize Data for JAGS        
jag.data <- list(n=data.bf.NY10.inf.smry$bloodmeals, temp=data.bf.NY10.inf.smry$temp,
                 N=data.bf.NY10.inf.smry$bloodmeals_offered,
                 N.obs=length(data.bf.NY10.inf.smry$temp))
jag.data

##### Run JAGS
#NY10.inf.bf <- jags(data=jag.data, inits=inits, parameters.to.save=parameters, 
#                       model.file="bloodfeeding_unexp_priors.txt", n.thin=nt, n.chains=nc, 
#                       n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd())
NY10.inf.bf <- do.call(jags.parallel, list(data=jag.data, inits=inits, parameters.to.save=parameters, 
                                             model.file="bloodfeeding_unexp_priors.txt", 
                                             n.thin=nt, n.chains=nc, 
                                             n.burnin=nb, n.iter=ni, DIC=TRUE, working.directory=getwd(),
                                             jags.seed = 42))

NY10.inf.bf
mcmcplot(NY10.inf.bf)


chains.NY10.inf.bf <- MCMCchains(NY10.inf.bf, params=c("Tmin", "Tmax", "rmax", 
                                                             "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.NY10.inf.bf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.NY10.inf.bf <- apply(curves, 1, mean)
CI.NY10.inf.bf <- apply(curves, 1, quantile, c(0.025, 0.975))


plot(temperatures, data.bf.NY10.inf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4), xlim=c(0, 40),
     main="NY10 (infected)")
arrows(temperatures, data.bf.NY10.inf.smry$prop, y1= data.bf.NY10.inf.smry$prop + data.bf.NY10.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.NY10.inf.smry$prop, y1= data.bf.NY10.inf.smry$prop - data.bf.NY10.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.inf.bf, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.inf.bf[1, ],
                                rev(CI.NY10.inf.bf[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)


unexp.bf.out
WN02.uninf.bf
WN02.inf.bf
NY10.uninf.bf
NY10.inf.bf

# Save MCMC output.
saveRDS(unexp.bf.out, './a_unexposed.RDS')
saveRDS(WN02.uninf.bf, './a_WN02_uninfected.RDS')
saveRDS(WN02.inf.bf, './a_WN02_infected.RDS')
saveRDS(NY10.uninf.bf, './a_NY10_uninfected.RDS')
saveRDS(NY10.inf.bf, './a_NY10_infected.RDS')


#View(unexp.bf.out$BUGSoutput$summary[c("Tmin", "Topt", "Tmax", "rmax", "alpha", "beta"), 
#                                      c(1,2,3,7,8,9)])
#View(WN02.uninf.bf$BUGSoutput$summary[c("Tmin", "Topt", "Tmax", "rmax", "alpha", "beta"), 
#c(1,2,3,7,8,9)])
#View(WN02.inf.bf$BUGSoutput$summary[c("Tmin", "Topt", "Tmax", "rmax", "alpha", "beta"), 
#                                      c(1,2,3,7,8,9)])

#View(NY10.uninf.bf$BUGSoutput$summary[c("Tmin", "Topt", "Tmax", "rmax", "alpha", "beta"), 
#                                      c(1,2,3,7,8,9)])
#View(NY10.inf.bf$BUGSoutput$summary[c("Tmin", "Topt", "Tmax", "rmax", "alpha", "beta"), 
#                                    c(1,2,3,7,8,9)])

par(mfrow=c(1, 1))
plot("", xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4), xlim=c(5, 35),
     main="Unexposed vs infected")
lines(temps, meancurve.unexp.bf, col="steelblue", lwd=2)
polygon(c(temps, rev(temps)), c(CI.unexp.bf[1, ],
                                rev(CI.unexp.bf[2,])), 
        col=alpha("steelblue", 0.3), lty=0)
lines(temps, meancurve.WN02.inf.bf, col="darkgreen", lwd=2)
polygon(c(temps, rev(temps)), c(CI.WN02.inf.bf[1, ],
                                rev(CI.WN02.inf.bf[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)
lines(temps, meancurve.NY10.inf.bf, col="darkblue", lwd=2)
polygon(c(temps, rev(temps)), c(CI.NY10.inf.bf[1, ],
                                rev(CI.NY10.inf.bf[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)


unexp.bf.out
WN02.inf.bf
NY10.inf.bf
