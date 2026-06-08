# Run this script from the repository root directory.

library('scales')
library('R2jags')
library('mcmcplots')
library('MCMCvis')

# Global parameters
errbar.length = 0.05
errbar.width = 1

WN02.name <- "2003.2"
NY10.name <- "2017.1"

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

# Overlay cycling proportion data (not used for fitting)
add_cycle_prop <- function(df, xcol = "temperature", ycol = "prop", secol = "sem",
                           pch = 17, col = "firebrick") {
  if(nrow(df) == 0) return(invisible(NULL))
  points(df[[xcol]], df[[ycol]], pch = pch, col = col)
  arrows(df[[xcol]], df[[ycol]], y1 = df[[ycol]] + df[[secol]],
         angle = 90, length = errbar.length, lwd = errbar.width, col = col)
  arrows(df[[xcol]], df[[ycol]], y1 = df[[ycol]] - df[[secol]],
         angle = 90, length = errbar.length, lwd = errbar.width, col = col)
}

# Overlay cycling continuous data (not used for fitting)
add_cycle_cont <- function(temp, mean, se, pch = 17, col = "firebrick") {
  keep <- is.finite(temp) & is.finite(mean) & is.finite(se)
  if(!any(keep)) return(invisible(NULL))
  points(temp[keep], mean[keep], pch = pch, col = col)
  arrows(temp[keep], mean[keep], y1 = mean[keep] + se[keep],
         angle = 90, length = errbar.length, lwd = errbar.width, col = col)
  arrows(temp[keep], mean[keep], y1 = mean[keep] - se[keep],
         angle = 90, length = errbar.length, lwd = errbar.width, col = col)
}



## Biting rate
data.bf <- read.csv("./data/bloodfeeding_lf.csv")
data.bf$infected <- as.logical(data.bf$infected)
data.bf$disseminated <- as.logical(data.bf$disseminated)
data.bf$cycle <- as.logical(data.bf$cycle)

data.bf$infected
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


# Infected biting rate, constant temperature
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

# Cycling BF data (not used for fitting)
data.bf.unexp.cyc <- subset(data.bf, ((data.bf$strain == "None") &
                                        (data.bf$cycle == TRUE) &
                                        (data.bf$infected == FALSE)))
data.bf.WN02.uninf.cyc <- subset(data.bf, ((data.bf$strain == "WN02") &
                                             (data.bf$cycle == TRUE) &
                                             (data.bf$infected == FALSE)))
data.bf.NY10.uninf.cyc <- subset(data.bf, ((data.bf$strain == "NY10") &
                                             (data.bf$cycle == TRUE) &
                                             (data.bf$infected == FALSE)))
data.bf.WN02.inf.cyc <- subset(data.bf, ((data.bf$strain == "WN02") &
                                           (data.bf$cycle == TRUE) &
                                           (data.bf$infected == TRUE)))
data.bf.NY10.inf.cyc <- subset(data.bf, ((data.bf$strain == "NY10") &
                                           (data.bf$cycle == TRUE) &
                                           (data.bf$infected == TRUE)))

summarize_bf_cycle <- function(df) {
  if(nrow(df) == 0) {
    return(data.frame(temperature=numeric(0), bloodmeals=numeric(0),
                      bloodmeals_offered=numeric(0), prop=numeric(0), sem=numeric(0)))
  }
  out <- aggregate(df[, c("bloodmeals", "bloodmeals_offered")],
                   by = list(temperature = df$temperature), FUN = sum)
  out$prop <- out$bloodmeals / out$bloodmeals_offered
  out$sem <- sqrt(out$prop * (1 - out$prop) / out$bloodmeals_offered)
  out
}

data.bf.unexp.cyc.smry <- summarize_bf_cycle(data.bf.unexp.cyc)
data.bf.WN02.uninf.cyc.smry <- summarize_bf_cycle(data.bf.WN02.uninf.cyc)
data.bf.NY10.uninf.cyc.smry <- summarize_bf_cycle(data.bf.NY10.uninf.cyc)
data.bf.WN02.inf.cyc.smry <- summarize_bf_cycle(data.bf.WN02.inf.cyc)
data.bf.NY10.inf.cyc.smry <- summarize_bf_cycle(data.bf.NY10.inf.cyc)


unexp.bf.out <- readRDS('./a_unexposed.RDS')
WN02.uninf.bf <- readRDS('./a_WN02_uninfected.RDS')
WN02.inf.bf <- readRDS('./a_WN02_infected.RDS')
NY10.uninf.bf <- readRDS('./a_NY10_uninfected.RDS')
NY10.inf.bf <- readRDS('./a_NY10_infected.RDS')


chains.unexp.bf <- MCMCchains(unexp.bf.out, params=c("Tmin", "Tmax", "rmax", 
                                                     "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.unexp.bf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.unexp.bf <- apply(curves, 1, mean)
CI.unexp.bf <- apply(curves, 1, quantile, c(0.025, 0.975))

par(mfrow=c(2, 3))

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
     main="2003.2 (exp. not infected)")
arrows(temperatures, data.bf.WN02.uninf.smry$prop, y1= data.bf.WN02.uninf.smry$prop + data.bf.WN02.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.WN02.uninf.smry$prop, y1= data.bf.WN02.uninf.smry$prop - data.bf.WN02.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)


# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.uninf.bf, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.uninf.bf[1, ],
                                rev(CI.WN02.uninf.bf[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)

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
     main="2017.1 (exp. not infected)")
arrows(temperatures, data.bf.NY10.uninf.smry$prop, y1= data.bf.NY10.uninf.smry$prop + data.bf.NY10.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.NY10.uninf.smry$prop, y1= data.bf.NY10.uninf.smry$prop - data.bf.NY10.uninf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)


# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.uninf.bf, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.uninf.bf[1, ],
                                rev(CI.NY10.uninf.bf[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)

chains.WN02.inf.bf <- MCMCchains(WN02.inf.bf, params=c("Tmin", "Tmax", "rmax", 
                                                       "alpha", "beta"))
temps <- seq(0, 40, 0.1)
curves <- apply(chains.WN02.inf.bf, 1,
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], x[5]))

# Find mean curve and credible intervals.
meancurve.WN02.inf.bf <- apply(curves, 1, mean)
CI.WN02.inf.bf <- apply(curves, 1, quantile, c(0.025, 0.975))

plot("")
plot(temperatures, data.bf.WN02.inf.smry$prop, pch=20, xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4), xlim=c(0, 40),
     main="2003.2 (infected)")
arrows(temperatures, data.bf.WN02.inf.smry$prop, y1= data.bf.WN02.inf.smry$prop + data.bf.WN02.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.WN02.inf.smry$prop, y1= data.bf.WN02.inf.smry$prop - data.bf.WN02.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.inf.bf, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.inf.bf[1, ],
                                rev(CI.WN02.inf.bf[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)

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
     main="2017.1 (infected)")
arrows(temperatures, data.bf.NY10.inf.smry$prop, y1= data.bf.NY10.inf.smry$prop + data.bf.NY10.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, data.bf.NY10.inf.smry$prop, y1= data.bf.NY10.inf.smry$prop - data.bf.NY10.inf.smry$sem, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.inf.bf, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.inf.bf[1, ],
                                rev(CI.NY10.inf.bf[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)


par(mfrow=c(1, 1))
plot("", xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4), xlim=c(5, 35),
     main="Biting rate", xaxt="n")
axis(1, at=seq(10, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 5:40, ticksize = -0.01, side = 1)
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


get_param_estimates <- function(x) {
  return(x[c("Tmin","Topt","Tmax","rmax","alpha","beta"), c(1,3,7)])
}
get_param_estimates(NY10.inf.bf$BUGSoutput$summary)

par(mfrow=c(1, 2))
plot("", xlim=c(0, 40), ylim=c(0, 6), xlab="Temperature [°C]", ylab="",
     main="Cardinal and optimal temperatures",
     yaxt = "n", xaxt="n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=1)
smry = list(unexp.bf.out$BUGSoutput$summary,
            WN02.uninf.bf$BUGSoutput$summary,
            WN02.inf.bf$BUGSoutput$summary,
            NY10.uninf.bf$BUGSoutput$summary,
            NY10.inf.bf$BUGSoutput$summary)
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

plot("", xlim=c(0, 0.5), ylim=c(0, 6), ylab="", xlab="bites / (mosquito * week)",
     main="Peak biting rate",
     yaxt = "n")

axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[4, 1], i, pch=20, col=colors[i])
  lines(c(params[4, 2], params[4, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

plot("", xlim=c(0, 1), ylim=c(0, 6), ylab="", xlab="alpha", yaxt = "n",
     main='Skewness')
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[5, 1], i, pch=20, col=colors[i])
  lines(c(params[5, 2], params[5, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

plot("", xlim=c(0, 1), ylim=c(0, 6), ylab="", xlab="beta", yaxt = "n",
     main='Breadth / tolerance ratio')
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[6, 1], i, pch=20, col=colors[i])
  lines(c(params[6, 2], params[6, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
}


par(mfrow=c(1, 2))
plot("", xlim=c(0, 40), ylim=c(0, 4), xlab="Temperature [°C]", ylab="",
     main="Cardinal and optimal temperatures",
     yaxt = "n", xaxt="n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)
axis(2, at=1:3, labels=c("unexp",  "2003.2",  "2017.1"), 
     las=1)
smry = list(unexp.bf.out$BUGSoutput$summary,
            WN02.inf.bf$BUGSoutput$summary,
            NY10.inf.bf$BUGSoutput$summary)
colors <- c("steelblue", "darkgreen", "darkblue")
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
text(10, 4 - eps, "Tmin")
text(24, 4 - eps, "Topt")
text(34, 4 - eps, "Tmax")

plot("", xlim=c(0, 0.5), ylim=c(0, 4), ylab="", xlab="bites / (mosquito * week)",
     main="Peak biting rate",
     yaxt = "n")

axis(2, at=1:3, labels=c("unexp",  "2003.2",  "2017.1"), 
     las=1)

colors <- c("steelblue", "darkgreen", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[4, 1], i, pch=20, col=colors[i])
  lines(c(params[4, 2], params[4, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}


## Lifespan

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

# Infected biting rate, constant temperature
data.lf.WN02.inf <- subset(data.lf, ((data.lf$genotype == "WN02") &
                                       (data.lf$cycle == FALSE) &
                                       (data.lf$infected == TRUE)))
data.lf.NY10.inf <- subset(data.lf, ((data.lf$genotype == "NY10") &
                                       (data.lf$cycle == FALSE) &
                                       (data.lf$infected == TRUE)))

unexp.lf.out.nb <- readRDS("./lf_unexposed.RDS")
WN02.uninf.lf.nb <- readRDS("./lf_WN02_uninfected.RDS")
NY10.uninf.lf.nb <- readRDS("./lf_NY10_uninfected.RDS")
WN02.inf.lf.nb <- readRDS("./lf_WN02_infected.RDS")
NY10.inf.lf.nb <- readRDS("./lf_NY10_infected.RDS")

chains.unexp.lf <- MCMCchains(unexp.lf.out.nb, params=c("Tmin", "Tmax", "rmax", 
                                                        "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.unexp.lf, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.unexp.lf <- apply(curves, 1, mean)
CI.unexp.lf <- apply(curves, 1, quantile, c(0.025, 0.975))

par(mfrow=c(2, 3))

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

chains.WN02.uninf.lf.nb <- MCMCchains(WN02.uninf.lf.nb, params=c("Tmin", "Tmax", "rmax", 
                                                                 "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.WN02.uninf.lf.nb, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.WN02.uninf.lf.nb <- apply(curves, 1, mean)
CI.WN02.uninf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))


m <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, mean))
s <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, sd))
n <- with(data.lf.WN02.uninf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="2003.2 (exp. not infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.uninf.lf.nb, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.uninf.lf.nb[1, ],
                                rev(CI.WN02.uninf.lf.nb[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)


chains.NY10.uninf.lf.nb <- MCMCchains(NY10.uninf.lf.nb, params=c("Tmin", "Tmax", "rmax", 
                                                                 "alpha", "beta"))


temps <- seq(0, 40, 0.1) 
curves <- apply(chains.NY10.uninf.lf.nb, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.NY10.uninf.lf.nb <- apply(curves, 1, mean)
CI.NY10.uninf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))


m <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, mean))
s <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, sd))
n <- with(data.lf.NY10.uninf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="2017.1 (exp. not infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.uninf.lf.nb, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.uninf.lf.nb[1, ],
                                rev(CI.NY10.uninf.lf.nb[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)


chains.WN02.inf.lf.nb <- MCMCchains(WN02.inf.lf.nb, params=c("Tmin", "Tmax", "rmax", 
                                                             "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.WN02.inf.lf.nb, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))

# Find mean curve and credible intervals.
meancurve.WN02.inf.lf.nb <- apply(curves, 1, mean)
CI.WN02.inf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))

m <- with(data.lf.WN02.inf, tapply(lifespan, temperature, mean))
s <- with(data.lf.WN02.inf, tapply(lifespan, temperature, sd))
n <- with(data.lf.WN02.inf, tapply(lifespan, temperature, length))

plot("")
plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="2003.2 (infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.inf.lf.nb, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.inf.lf.nb[1, ],
                                rev(CI.WN02.inf.lf.nb[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)


chains.NY10.inf.lf.nb <- MCMCchains(NY10.inf.lf.nb, params=c("Tmin", "Tmax", "rmax", 
                                                             "alpha", "beta"))

temps <- seq(0, 40, 0.1) 
curves <- apply(chains.NY10.inf.lf.nb, 1, 
                function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                    x[5]))
# Find mean curve and credible intervals.
meancurve.NY10.inf.lf.nb <- apply(curves, 1, mean)
CI.NY10.inf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))

m <- with(data.lf.NY10.inf, tapply(lifespan, temperature, mean))
s <- with(data.lf.NY10.inf, tapply(lifespan, temperature, sd))
n <- with(data.lf.NY10.inf, tapply(lifespan, temperature, length))

plot(temperatures, m,
     col=alpha("black", 1.0), xlim=c(0, 40), ylim=c(0, 40),
     xlab="Temperature [°C]", ylab="Lifespan [days]",
     main="2017.1 (infected)", pch=20)
arrows(temperatures, m, y1= m + s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(temperatures, m,  y1=m - s / sqrt(n), 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.inf.lf.nb, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.inf.lf.nb[1, ],
                                rev(CI.NY10.inf.lf.nb[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)

### Plot probabilities

curves.unexp <- apply(chains.unexp.lf, 1, 
                      function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                          x[5]))


curves.WN02 <- apply(chains.WN02.inf.lf.nb, 1, 
                     function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                         x[5]))

curves.NY10 <- apply(chains.NY10.inf.lf.nb, 1, 
                     function(x) flexTPC(temps, x[1], x[2], x[3], x[4], 
                                         x[5]))

par(mfrow=c(1,1))
plot("",
     col=alpha("black", 1.0), xlim=c(12, 35), ylim=c(-10, 10),
     xlab="Temperature [°C]", ylab="Lifespan infected - unexposed [days]",
     main="Lifespan difference", pch=20)

# Find mean curve and credible intervals.
curves <- curves.WN02 - curves.unexp
meancurve <- apply(curves, 1, mean)
CI <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot mean curve and 95% credible interval.
lines(temps, meancurve, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI[1, ],
                                rev(CI[2, ])), 
        col=alpha("darkgreen", 0.3), lty=0)

# Find mean curve and credible intervals.
curves <- curves.NY10 - curves.unexp
meancurve <- apply(curves, 1, mean)
CI <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot mean curve and 95% credible interval.
lines(temps, meancurve, col="darkblue")
polygon(c(temps, rev(temps)), c(CI[1, ],
                                rev(CI[2, ])), 
        col=alpha("darkblue", 0.3), lty=0)

abline(h=0, lty=2)

par(mfrow=c(1,1))
plot("",
     col=alpha("black", 1.0), xlim=c(12, 35), ylim=c(0.8, 1),
     xlab="Temperature [°C]", ylab="probability",
     main="Pr(Lf_inf < Lf_unexp)", pch=20)

# Find mean curve and credible intervals.
curves <- curves.WN02 < curves.unexp    
meancurve <- apply(curves, 1, mean)
#CI <- apply(curves, 1, quantile, c(0.025, 0.975), na.rm=TRUE)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve, col="darkgreen", lwd=2)
#polygon(c(temps, rev(temps)), c(CI[1, ],
#                                rev(CI[2, ])), 
#        col=alpha("darkgreen", 0.3), lty=0)

# Find mean curve and credible intervals.
curves <- curves.NY10 < curves.unexp
meancurve <- apply(curves, 1, mean)
#CI <- apply(curves, 1, quantile, c(0.025, 0.975))

# Plot mean curve and 95% credible interval.
lines(temps, meancurve, col="darkblue", lwd=2)


abline(h=0.95, lty=2)
abline(h=0.99, lty=2)


# Find mean curve and credible intervals.
meancurve.NY10.uninf.lf.nb <- apply(curves, 1, mean)
CI.NY10.uninf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))

# Find mean curve and credible intervals.
meancurve.NY10.uninf.lf.nb <- apply(curves, 1, mean)
CI.NY10.uninf.lf.nb <- apply(curves, 1, quantile, c(0.025, 0.975))


par(mfrow=c(1,1))

plot(temps, meancurve.unexp.lf, col="steelblue", type='l', xlim=c(0, 40), 
     ylim=c(0, 30),
     xlab="Temperature [°C]", ylab="Lifespan [days]", lwd=2,
     main="Lifespan", xaxt="n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(0, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)
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

par(mfrow=c(1, 2))
plot("", xlim=c(0, 40), ylim=c(0, 6), xlab="Temperature [°C]", ylab="",
     main="Cardinal and optimal temperatures", xaxt="n",
     yaxt = "n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)

axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=1)
smry = list(unexp.lf.out.nb$BUGSoutput$summary,
            WN02.uninf.lf.nb$BUGSoutput$summary,
            NY10.uninf.lf.nb$BUGSoutput$summary,
            WN02.inf.lf.nb$BUGSoutput$summary,
            NY10.inf.lf.nb$BUGSoutput$summary)
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
text(5, 6 - eps, "Tmin")
text(18, 6 - eps, "Topt")
text(34, 6 - eps, "Tmax")

plot("", xlim=c(0, 30), ylim=c(0, 6), ylab="", xlab="Lifespan [days]",
     main="Peak lifespan",
     yaxt = "n")
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[4, 1], i, pch=20, col=colors[i])
  lines(c(params[4, 2], params[4, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

plot("", xlim=c(0, 1), ylim=c(0, 6), ylab="", xlab="alpha", yaxt = "n",
     main='Skewness')
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[5, 1], i, pch=20, col=colors[i])
  lines(c(params[5, 2], params[5, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

plot("", xlim=c(0, 1), ylim=c(0, 6), ylab="", xlab="beta", yaxt = "n",
     main='Breadth / tolerance ratio')
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[6, 1], i, pch=20, col=colors[i])
  lines(c(params[6, 2], params[6, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}


par(mfrow=c(1, 2))
plot("", xlim=c(0, 40), ylim=c(0, 4), xlab="Temperature [°C]", ylab="",
     main="Cardinal and optimal temperatures", xaxt="n",
     yaxt = "n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)

axis(2, at=1:3, labels=c("unexp", "2003.2", "2017.1"), 
     las=1)
smry = list(unexp.lf.out.nb$BUGSoutput$summary,
            WN02.inf.lf.nb$BUGSoutput$summary,
            NY10.inf.lf.nb$BUGSoutput$summary)
colors <- c("steelblue", "darkgreen", "darkblue")
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
text(5, 4 - eps, "Tmin")
text(18, 4 - eps, "Topt")
text(34, 4 - eps, "Tmax")

plot("", xlim=c(0, 30), ylim=c(0, 4), ylab="", xlab="Lifespan [days]",
     main="Peak lifespan",
     yaxt = "n")
axis(2, at=1:3, labels=c("unexp", "2003.2", "2017.1"), 
     las=2)

colors <- c("steelblue", "darkgreen",  "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[4, 1], i, pch=20, col=colors[i])
  lines(c(params[4, 2], params[4, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}


## Oviposition

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

unexp.pO.out <- readRDS("./pO_unexposed.RDS")
WN02.uninf.pO.out <- readRDS("./pO_WN02_uninfected.RDS")
NY10.uninf.pO.out <- readRDS("./pO_NY10_uninfected.RDS")
WN02.inf.pO.out <- readRDS("./pO_WN02_infected.RDS")
NY10.inf.pO.out <- readRDS("./pO_NY10_infected.RDS")

## Plot TPCs

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
     main="2003.2 (exp. not infected)")
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
     main="2003.2 (infected)")
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
     main="Oviposition", xaxt="n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)
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


par(mfrow=c(1,1))

plot(temps, meancurve.pO.unexp, col="steelblue", type='l', xlim=c(0, 40), 
     ylim=c(0, 1),
     xlab="Temperature [°C]", ylab="prop. ovipositing", lwd=2,
     main="Oviposition", xaxt="n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)

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



par(mfrow=c(1, 2))
plot("", xlim=c(0, 40), ylim=c(0, 6), xlab="Temperature [°C]", ylab="",
     main="Cardinal and optimal temperatures", xaxt="n",
     yaxt = "n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
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
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[4, 1], i, pch=20, col=colors[i])
  lines(c(params[4, 2], params[4, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

plot("", xlim=c(0, 1), ylim=c(0, 6), ylab="", xlab="alpha", yaxt = "n",
     main='Skewness')
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[5, 1], i, pch=20, col=colors[i])
  lines(c(params[5, 2], params[5, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

plot("", xlim=c(0, 1), ylim=c(0, 6), ylab="", xlab="beta", yaxt = "n",
     main='Breadth / tolerance ratio')
axis(2, at=1:5, labels=c("unexp", "2003.2n", "2003.2i", "2017.1n", "2017.1i"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkgreen", "darkblue", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[6, 1], i, pch=20, col=colors[i])
  lines(c(params[6, 2], params[6, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

par(mfrow=c(1, 2))
plot("", xlim=c(0, 40), ylim=c(0, 4), xlab="Temperature [°C]", ylab="",
     main="Cardinal and optimal temperatures", xaxt="n",
     yaxt = "n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)
axis(2, at=1:3, labels=c("unexp", "2003.2", "2017.1"), 
     las=1)
smry = list(unexp.pO.out$BUGSoutput$summary,
            WN02.inf.pO.out$BUGSoutput$summary,
            NY10.inf.pO.out$BUGSoutput$summary)
colors <- c("steelblue", "darkgreen", "darkblue")
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
text(10, 4 - eps, "Tmin")
text(24, 4 - eps, "Topt")
text(34, 4 - eps, "Tmax")

plot("", xlim=c(0, 1), ylim=c(0, 4), ylab="", xlab="prob. ovipositing",
     main="Peak oviposition",
     yaxt = "n")
axis(2, at=1:3, labels=c("unexp", "2003.2", "2017.1"), 
     las=2)

colors <- c("steelblue", "darkgreen", "darkblue")
for(i in 1:length(smry)) {
  params <- get_param_estimates(smry[[i]])
  #print(params)
  points(params[4, 1], i, pch=20, col=colors[i])
  lines(c(params[4, 2], params[4, 3]), rep(i, 2), lwd=2,
        col=alpha(colors[i], 0.5))
  
}

###### Vector competence.

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

#errbar.length=0.06
#errbar.width=2
temps <- seq(0, 40, 0.1)

mech.out <- readRDS("bbc_mech.RDS")
chains.WN02.inf <- MCMCchains(mech.out, params=c("im_T50", "im_n", "im_rmin", 
                                                 "vr_WN02_Tmin", "vr_WN02_Tmax", "vr_WN02_rmax", "vr_WN02_alpha", "vr_WN02_beta",
                                                 "b_WN02_v50", "b_WN02_n"))

curves <- apply(chains.WN02.inf, 1, 
                function(x) bbc(temps, x[["im_rmin"]], 1.0, x[["im_T50"]], x[["im_n"]], 
                                x[["vr_WN02_Tmin"]], x[["vr_WN02_Tmax"]], x[["vr_WN02_rmax"]], x[["vr_WN02_alpha"]], x[["vr_WN02_beta"]],
                                x[["b_WN02_n"]]))

# Find mean curve and credible intervals.
meancurve.WN02.inf <- apply(curves, 1, mean)
CI.WN02.inf <- apply(curves, 1, quantile, c(0.025, 0.975))

par(mfrow=c(2,2))

plot(data.WN02$temperature, data.WN02$p.infected, main="2003.2",
     xlab="Temperature [°C]", ylab="b (prop. infected)", xlim=c(0, 40),
     ylim=c(0,1), pch=20, xaxt="n")
axis(1, at=seq(0, 40, 10))
rug(x = seq(5, 35, 5), ticksize = -0.03, side = 1)
rug(x = 0:40, ticksize = -0.01, side = 1)

arrows(data.WN02$temperature, data.WN02$p.infected, y1=data.WN02$p.infected+data.WN02$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.WN02$temperature, data.WN02$p.infected, y1=data.WN02$p.infected-data.WN02$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.inf, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.inf[1,],
                                rev(CI.WN02.inf[2,])), 
        col=alpha("darkgreen", 0.3), lty=0)


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

plot(data.NY10$temperature, data.NY10$p.infected, main="2017.1",
     xlab="Temperature [°C]", ylab="b (prop. infected)", xlim=c(0, 40),
     ylim=c(0,1), pch=20)
arrows(data.NY10$temperature, data.NY10$p.infected, y1=data.NY10$p.infected+data.NY10$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.NY10$temperature, data.NY10$p.infected, y1=data.NY10$p.infected-data.NY10$p.infected.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.inf, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.inf[1,],
                                rev(CI.NY10.inf[2,])), 
        col=alpha("darkblue", 0.3), lty=0)


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

plot(data.WN02$temperature, data.WN02$p.disseminated, main="2003.2",
     xlab="Temperature [°C]", ylab="bc (prop. disseminated)", xlim=c(0, 40),
     ylim=c(0, 1), pch=17)
arrows(data.WN02$temperature, data.WN02$p.disseminated, y1=data.WN02$p.disseminated +data.WN02$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.WN02$temperature, data.WN02$p.disseminated, y1=data.WN02$p.disseminated-data.WN02$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)

# Plot mean curve and 95% credible interval.
lines(temps, meancurve.WN02.dis, col="darkgreen")
polygon(c(temps, rev(temps)), c(CI.WN02.dis[1,],
                                rev(CI.WN02.dis[2,])), 
        col=alpha("darkgreen", 0.3), lty=0)

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


plot(data.NY10$temperature, data.NY10$p.disseminated, main="2017.1",
     xlab="Temperature [°C]", ylab="bc (prop. disseminated)", xlim=c(0, 40),
     ylim=c(0,1), pch=17)
arrows(data.NY10$temperature, data.NY10$p.disseminated, y1=data.NY10$p.disseminated+data.NY10$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)
arrows(data.NY10$temperature, data.NY10$p.disseminated, y1=data.NY10$p.disseminated-data.NY10$p.disseminated.sderr, 
       angle=90, length=errbar.length, lwd=errbar.width)


# Plot mean curve and 95% credible interval.
lines(temps, meancurve.NY10.dis, col="darkblue")
polygon(c(temps, rev(temps)), c(CI.NY10.dis[1,],
                                rev(CI.NY10.dis[2,])), 
        col=alpha("darkblue", 0.3), lty=0)



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
legend(5, 4, c("immunity", "replication 2003.2", "replication 2017.1"), col=c("gold", "darkgreen", "darkblue"),
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
legend(6, 0.8, c("immunity", "replication 2003.2", "replication 2017.1"), col=c("gold", "darkgreen", "darkblue"),
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
     ylim=c(-0.1, 1.1), type='l', main="2003.2")
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
     ylim=c(-0.1, 1.1), type='l', main="2017.1")
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