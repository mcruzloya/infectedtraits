# Uses fitted trait models to calculate and plot relative R0 curves.
# Run from the repository root directory.

shocket_traits_dir <- "./shocket_fits/"

set.seed(42) # Set seed for reproducibility.

# We will need these packages. Please make sure you install them beforehand!
library('R2jags')
library('mcmcplots')
library('MCMCvis')
library('scales')

logM <- function(aG, lf, EFGC, EV, pLA, MDR) {
  # aG: Inverse gonotrophic cycle length.
  eps <- 1e-10 # Small constant to prevent denominators being zero.
  logM <- log(aG) + log(EFGC) + log(EV) + log(pLA) + log(MDR) + 2 * log(lf)
  return(logM)
}

# Relative R0 without specific parameters for infectious mosquitoes.
R0 <- function(a, aG, bc, lf, PDR, EFGC, EV, pLA, MDR) {
  eps <- 1e-10 # Small constant to prevent denominators being zero 
  # Calculate as log of R0^2 for numerical stability.
  logR02 <- (2 * log(a) + log(bc) - (1 / (PDR * lf + eps)) + log(lf) 
             + logM(aG, lf, EFGC, EV, pLA, MDR))
  return(sqrt(exp(logR02)))
}

# Relative R0 with specific parameters for infectious mosquitoes.
R0_I <- function(a, aI, aG, bc, lf, lfI, PDR, EFGC, EV, pLA, MDR) {
  eps <- 1e-10 # Small constant to prevent denominators being zero 
  # Calculate as log of R0^2 for numerical stability.
  logR02 <- (log(a) + log(aI) + log(bc) - (1 / (PDR * lfI + eps)) + log(lfI) 
             + logM(aG, lf, EFGC, EV, pLA, MDR))
  return(sqrt(exp(logR02)))
}

# Briere TPC model.
briere <- function(T, Tmin, Tmax, c) {
  result <- c * T * (T - Tmin) * sqrt((Tmax - T) * (T > Tmin) * (T < Tmax))
  return(result)
}

# Quadratic TPC model.
quad <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(result)
}

# Quadratic TPC model with upper limit at 1 for proportions.
quad_lim <- function(T, Tmin, Tmax, c) {
  result <- (T > Tmin) * (T < Tmax) * c * (T - Tmin) * (Tmax - T)
  return(pmin(result, 1))
}

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

plot_curves <- function(temps, curves, xlim=c(0, 50), ylim=c(0, 1), col='steelblue',
                        a=0.3, lwd=2) {
  meancurve <- apply(curves, 1, mean)
  CI <- apply(curves, 1, quantile, c(0.025, 0.975))
  lines(temps, meancurve, type="l", col=col, xlim=xlim, ylim=ylim, lwd=2)
  polygon(c(temps, rev(temps)), c(CI[1, ],
                                  rev(CI[2, ])), 
          col=alpha(col, a), lty=0)
}

temps <- seq(0, 50, 0.05) # Test with small vector first to debug.

# Indices of MCMC chains to use for generating the R0 curves.
idx <- seq(1, 250000, 5)

## Load data for traits in unexposed/exposed experiments.

# Biting rate
a.unexp <- readRDS("./a_unexposed.RDS")
a.WN02.inf <- readRDS("./a_WN02_infected.RDS")
a.NY10.inf <- readRDS("./a_NY10_infected.RDS")

a.unexp.chains <- MCMCchains(a.unexp, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))[idx, ]
a.WN02.inf.chains <- MCMCchains(a.WN02.inf, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))[idx, ]
a.NY10.inf.chains <- MCMCchains(a.NY10.inf, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))[idx, ]

a.unexp.curves <- apply(a.unexp.chains, 1, 
                        function(x) flexTPC(temps, x[["Tmin"]], x[["Tmax"]],
                                            x[["rmax"]], x[["alpha"]],
                                            x[["beta"]]))
a.WN02.inf.curves <- apply(a.WN02.inf.chains, 1, 
                           function(x) flexTPC(temps, x[["Tmin"]], x[["Tmax"]],
                                               x[["rmax"]], x[["alpha"]],
                                               x[["beta"]]))
a.NY10.inf.curves <- apply(a.NY10.inf.chains, 1, 
                           function(x) flexTPC(temps, x[["Tmin"]], x[["Tmax"]],
                                               x[["rmax"]], x[["alpha"]],
                                               x[["beta"]]))
par(mfrow=c(1, 1))
plot("", xlab='Temperature [°C]',
     ylab="Biting rate [bites / (mosquito * week)]", ylim=c(0, 0.4), xlim=c(0, 40),
     main="Biting rate")
plot_curves(temps, a.unexp.curves, col="steelblue")
plot_curves(temps, a.WN02.inf.curves, col="darkgreen")
plot_curves(temps, a.NY10.inf.curves, col="darkblue")



# Vector competence.
b_bc <- readRDS("./bbc_mech.RDS")
b_bc
bc.chains <- MCMCchains(b_bc, params=c("im_T50", "im_n", "im_rmin", 
                                       "vr_WN02_Tmin", "vr_WN02_Tmax", "vr_WN02_rmax", #"vr_WN02_alpha", "vr_WN02_beta",
                                       "bc_WN02_v50", "bc_WN02_n",
                                       "vr_NY10_Tmin", "vr_NY10_Tmax", "vr_NY10_rmax", # "vr_NY10_alpha", "vr_NY10_beta",
                                       "bc_NY10_v50", "bc_NY10_n"))[idx, ]
bc.WN02.curves <- apply(bc.chains, 1, 
                        function(x) bbc(temps, x[["im_rmin"]], 1.0, x[["im_T50"]], x[["im_n"]], 
                                        x[["vr_WN02_Tmin"]], x[["vr_WN02_Tmax"]], x[["vr_WN02_rmax"]], 0.8, 0.2,
                                        x[["bc_WN02_n"]]))
bc.NY10.curves <- apply(bc.chains, 1, 
                        function(x) bbc(temps, x[["im_rmin"]], 1.0, x[["im_T50"]], x[["im_n"]], 
                                        x[["vr_NY10_Tmin"]], x[["vr_NY10_Tmax"]], x[["vr_NY10_rmax"]], 0.8, 0.2,
                                        x[["bc_NY10_n"]]))

# Lifespan
lf.unexp <- readRDS("./lf_unexposed.RDS")
lf.WN02.inf <- readRDS("./lf_WN02_infected.RDS")
lf.NY10.inf <- readRDS("./lf_NY10_infected.RDS")


# Oviposition
pO.unexp <- readRDS("./pO_unexposed.RDS")
pO.WN02.inf <- readRDS("./pO_WN02_infected.RDS")
pO.NY10.inf <- readRDS("./pO_NY10_infected.RDS")


plot_temps <- function(smry_list, colors=c("steelblue", "darkgreen", "darkblue"),
                       ylab=c("U", "WN02", "NY10"), main="Lifespan") {
  plot("", xlim=c(0, 50), ylim=c(0.5, length(smry_list) + 0.5), xlab="Temperature [°C]", ylab="", main=main,
       yaxt="n")
  axis(2, at = 1:length(smry_list), labels=ylab, las=2)
  for(i in 1:length(smry_list)) {
    smry = smry_list[[i]]
    col = colors[i]
    
    for(param in c("Tmin", "Topt", "Tmax")){
      m <- smry[param, "mean"]
      CI <- c(smry[param, "2.5%"], smry[param, "97.5%"])
      points(m, i, col=col, pch=20)
      lines(CI, c(i, i), col=col, lwd=2)
    }
  }
}


plot_temps(list(lf.unexp$BUGSoutput$summary, lf.WN02.inf$BUGSoutput$summary, lf.NY10.inf$BUGSoutput$summary),
           main="Lifespan")
plot_temps(list(a.unexp$BUGSoutput$summary, a.WN02.inf$BUGSoutput$summary, a.NY10.inf$BUGSoutput$summary),
           main="Biting rate")
plot_temps(list(pO.unexp$BUGSoutput$summary, pO.WN02.inf$BUGSoutput$summary, pO.NY10.inf$BUGSoutput$summary),
           main="Oviposition")


lf.unexp.chains <- MCMCchains(lf.unexp, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))[idx, ]
lf.WN02.inf.chains <- MCMCchains(lf.WN02.inf, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))[idx, ]
lf.NY10.inf.chains <- MCMCchains(lf.NY10.inf, params=c("Tmin", "Tmax", "rmax", "alpha", "beta"))[idx, ]



lf.unexp.curves <- apply(lf.unexp.chains, 1, 
                         function(x) flexTPC(temps, x[["Tmin"]], x[["Tmax"]],
                                             x[["rmax"]], x[["alpha"]],
                                             x[["beta"]]))
lf.WN02.inf.curves <- apply(lf.WN02.inf.chains, 1, 
                            function(x) flexTPC(temps, x[["Tmin"]], x[["Tmax"]],
                                                x[["rmax"]], x[["alpha"]],
                                                x[["beta"]]))
lf.NY10.inf.curves <- apply(lf.NY10.inf.chains, 1, 
                            function(x) flexTPC(temps, x[["Tmin"]], x[["Tmax"]],
                                                x[["rmax"]], x[["alpha"]],
                                                x[["beta"]]))

par(mfrow=c(1, 1))
plot("", xlab='Temperature [°C]',
     ylab="Lifespan [days]", ylim=c(0, 30), xlim=c(0, 40),
     main="Lifespan")
plot_curves(temps, lf.unexp.curves, col="steelblue")
plot_curves(temps, lf.WN02.inf.curves, col="darkgreen")
plot_curves(temps, lf.NY10.inf.curves, col="darkblue")


### Additional trait fits used in R0 calculations.

idx_trait <- round(seq(1, 80000, 80000/ 50000))
idx_trait
load(paste(shocket_traits_dir, "jagsout_a_Cpip_inf.Rdata", sep=""))
aG.chains <- MCMCchains(a.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))[]
aG.chains

aG.curves <-  apply(aG.chains[idx_trait, ], 1, function(x) briere(temps, x[1], x[2], x[3]))

load(paste(shocket_traits_dir, "jagsout_EFOC_Cpip_inf.Rdata", sep=""))
EFOC.chains <- MCMCchains(EFOC.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
EFOC.curves <- apply(EFOC.chains[idx_trait, ], 1, function(x) quad(temps, x[1], x[2], x[3]))

load(paste(shocket_traits_dir, "jagsout_EV_Cpip_inf.Rdata", sep=""))
EV.chains <- MCMCchains(EV.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
EV.curves <- apply(EV.chains[idx_trait, ], 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

load(paste(shocket_traits_dir, "jagsout_MDR_Cpip_inf.Rdata", sep=""))
MDR.chains <- MCMCchains(MDR.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
MDR.curves <- apply(MDR.chains[idx_trait, ], 1, function(x) briere(temps, x[1], x[2], x[3]))

load(paste(shocket_traits_dir, "jagsout_PDR_CpipWNV_inf.Rdata", sep=""))
PDR.chains <- MCMCchains(PDR.CpipWNV.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
PDR.curves <- apply(PDR.chains[idx_trait, ], 1, function(x) briere(temps, x[1], x[2], x[3]))

load(paste(shocket_traits_dir, "jagsout_pLA_Cpip_inf.Rdata", sep=""))
pLA.chains <- MCMCchains(pLA.Cpip.out.inf, params=c("cf.T0", "cf.Tm", "cf.q"))
pLA.curves <- apply(pLA.chains[idx_trait, ], 1, function(x) quad_lim(temps, x[1], x[2], x[3]))

str(pLA.curves)

# Thin bc curves to make same size as others.
#split_func <- function(x, by) {
#  r <- diff(range(x))
#  out <- seq(0, r - by - 1, by = by)
#  c(round(min(x) + c(0, out - 0.51 + (max(x) - max(out)) / 2), 0), max(x))
#}

#idx <- split_func(1:200000, 200000 / 75000)[2:75001]
#str(idx)
#str(bc.WN02.curves)
#bc.WN02.curves <- bc.WN02.curves[ , idx]
#bc.NY10.curves <- bc.NY10.curves[ , idx]

#str(bc.WN02.curves)

#### R0 calculations

find_Topt <- function(temps, curves) {
  idx <- apply(curves, 2, which.max)
  return(temps[idx])
}
find_Tmin <- function(temps, curves) {
  idx <- apply(curves, 2, function(x) which(x > 0)[1] - 1)
  return(temps[idx])
}

find_Tmax <- function(temps, curves) {
  idx <- apply(curves, 2, function(x) tail(which(x > 0), 1) + 1)
  return(temps[idx])
}

## R0 calculations

#WN02
R0.WN02.curves <- R0(a.unexp.curves / 7, aG.curves, bc.WN02.curves, lf.unexp.curves,
                     PDR.curves, EFOC.curves, EV.curves, pLA.curves, MDR.curves)
R0_I.WN02.curves <- R0_I(a.unexp.curves / 7, a.WN02.inf.curves / 7, aG.curves, bc.WN02.curves,
                         lf.unexp.curves, lf.WN02.inf.curves, PDR.curves, EFOC.curves,
                         EV.curves, pLA.curves, MDR.curves)

#NY10
R0.NY10.curves <- R0(a.unexp.curves / 7, aG.curves, bc.NY10.curves, lf.unexp.curves,
                     PDR.curves, EFOC.curves, EV.curves, pLA.curves, MDR.curves)
R0_I.NY10.curves <- R0_I(a.unexp.curves / 7, a.NY10.inf.curves / 7, aG.curves, bc.NY10.curves,
                         lf.unexp.curves, lf.NY10.inf.curves, PDR.curves, EFOC.curves,
                         EV.curves, pLA.curves, MDR.curves)



plot("", main="Relative R0", xlab="Temperature [°C]", ylab="Relative R0",
     xlim=c(10, 40), ylim=c(0, 6.0),
     xaxt="n", yaxt="n")
axis(1, at=seq(10, 40, 10))
rug(x = seq(15, 35, 5), ticksize = -0.03, side = 1)
rug(x = 10:40, ticksize = -0.01, side = 1)

meancurve.WN02.R0 <- apply(R0.WN02.curves, 1, mean)
CI.WN02.R0 <- apply(R0.WN02.curves, 1, quantile, c(0.025, 0.975))

meancurve.NY10.R0 <- apply(R0.NY10.curves, 1, mean)
CI.NY10.R0 <- apply(R0.NY10.curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.WN02.R0, col="darkgreen", lwd=2)
polygon(c(temps, rev(temps)), c(CI.WN02.R0[1,], rev(CI.WN02.R0[2,])), 
        col=alpha("darkgreen", 0.3), lty=0)

lines(temps, meancurve.NY10.R0, col="darkblue", lwd=2)
polygon(c(temps, rev(temps)), c(CI.NY10.R0[1,], rev(CI.NY10.R0[2,])), 
        col=alpha("darkblue", 0.3), lty=0)

plot("", main="Relative R0 (infected param)", xlab="Temperature [°C]", ylab="Relative R0",
     xlim=c(10, 40), ylim=c(0, 6.0),
     xaxt="n", yaxt="n")
axis(1, at=seq(10, 40, 10))
rug(x = seq(15, 35, 5), ticksize = -0.03, side = 1)
rug(x = 10:40, ticksize = -0.01, side = 1)

meancurve.WN02.R0_I <- apply(R0_I.WN02.curves, 1, mean)
CI.WN02.R0_I <- apply(R0_I.WN02.curves, 1, quantile, c(0.025, 0.975))

meancurve.NY10.R0_I <- apply(R0_I.NY10.curves, 1, mean)
CI.NY10.R0_I <- apply(R0_I.NY10.curves, 1, quantile, c(0.025, 0.975))

lines(temps, meancurve.WN02.R0_I, col="darkgreen", lwd=2)
polygon(c(temps, rev(temps)), c(CI.WN02.R0_I[1,], rev(CI.WN02.R0_I[2,])), 
        col=alpha("darkgreen", 0.3), lty=0)

lines(temps, meancurve.NY10.R0_I, col="darkblue", lwd=2)
polygon(c(temps, rev(temps)), c(CI.NY10.R0_I[1,], rev(CI.NY10.R0_I[2,])), 
        col=alpha("darkblue", 0.3), lty=0)

R0.summary <- data.frame(params = c(rep('unexp', 2), rep('inf', 2)),
                         strain=rep(c("WN02", "NY10"), 2),
                         Tmin=rep(0, 4), Tminl=rep(0,4), Tminh=rep(0,4),
                         Topt=rep(0, 4), Toptl=rep(0,4), Topth=rep(0,4),
                         Tmax=rep(0, 4), Tmaxl=rep(0,4), Tmaxh=rep(0,4)
)
R0.summary


## WN02 Unexposed
Topt.chains <- find_Topt(temps, R0.WN02.curves)
Topt.CI <- quantile(Topt.chains, c(0.025, 0.975))

R0.summary$Topt[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- mean(Topt.chains)
R0.summary$Toptl[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- Topt.CI[1]
R0.summary$Topth[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- Topt.CI[2]

Tmin.chains <- find_Tmin(temps, R0.WN02.curves)
Tmin.CI <- quantile(Tmin.chains, c(0.025, 0.975))

R0.summary$Tmin[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- mean(Tmin.chains)
R0.summary$Tminl[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- Tmin.CI[1]
R0.summary$Tminh[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- Tmin.CI[2]

Tmax.chains <- find_Tmax(temps, R0.WN02.curves)
Tmax.CI <- quantile(Tmax.chains, c(0.025, 0.975))

R0.summary$Tmax[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- mean(Tmax.chains)
R0.summary$Tmaxl[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- Tmax.CI[1]
R0.summary$Tmaxh[(R0.summary$strain == "WN02") & (R0.summary$params == "unexp")] <- Tmax.CI[2]

## NY10 Unexposed
Topt.chains <- find_Topt(temps, R0.NY10.curves)
Topt.CI <- quantile(Topt.chains, c(0.025, 0.975))

R0.summary$Topt[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- mean(Topt.chains)
R0.summary$Toptl[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- Topt.CI[1]
R0.summary$Topth[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- Topt.CI[2]

Tmin.chains <- find_Tmin(temps, R0.NY10.curves)
Tmin.CI <- quantile(Tmin.chains, c(0.025, 0.975))

R0.summary$Tmin[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- mean(Tmin.chains)
R0.summary$Tminl[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- Tmin.CI[1]
R0.summary$Tminh[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- Tmin.CI[2]

Tmax.chains <- find_Tmax(temps, R0.NY10.curves)
Tmax.CI <- quantile(Tmax.chains, c(0.025, 0.975))

R0.summary$Tmax[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- mean(Tmax.chains)
R0.summary$Tmaxl[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- Tmax.CI[1]
R0.summary$Tmaxh[(R0.summary$strain == "NY10") & (R0.summary$params == "unexp")] <- Tmax.CI[2]

## WN02 Infected
Topt.chains <- find_Topt(temps, R0_I.WN02.curves)
Topt.CI <- quantile(Topt.chains, c(0.025, 0.975))

R0.summary$Topt[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- mean(Topt.chains)
R0.summary$Toptl[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- Topt.CI[1]
R0.summary$Topth[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- Topt.CI[2]

Tmin.chains <- find_Tmin(temps, R0_I.WN02.curves)
Tmin.CI <- quantile(Tmin.chains, c(0.025, 0.975))

R0.summary$Tmin[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- mean(Tmin.chains)
R0.summary$Tminl[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- Tmin.CI[1]
R0.summary$Tminh[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- Tmin.CI[2]

Tmax.chains <- find_Tmax(temps, R0_I.WN02.curves)
Tmax.CI <- quantile(Tmax.chains, c(0.025, 0.975))

R0.summary$Tmax[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- mean(Tmax.chains)
R0.summary$Tmaxl[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- Tmax.CI[1]
R0.summary$Tmaxh[(R0.summary$strain == "WN02") & (R0.summary$params == "inf")] <- Tmax.CI[2]

## NY10 Infected
Topt.chains <- find_Topt(temps, R0_I.NY10.curves)
Topt.CI <- quantile(Topt.chains, c(0.025, 0.975))

R0.summary$Topt[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- mean(Topt.chains)
R0.summary$Toptl[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- Topt.CI[1]
R0.summary$Topth[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- Topt.CI[2]

Tmin.chains <- find_Tmin(temps, R0_I.NY10.curves)
Tmin.CI <- quantile(Tmin.chains, c(0.025, 0.975))

R0.summary$Tmin[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- mean(Tmin.chains)
R0.summary$Tminl[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- Tmin.CI[1]
R0.summary$Tminh[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- Tmin.CI[2]

Tmax.chains <- find_Tmax(temps, R0_I.NY10.curves)
Tmax.CI <- quantile(Tmax.chains, c(0.025, 0.975))

R0.summary$Tmax[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- mean(Tmax.chains)
R0.summary$Tmaxl[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- Tmax.CI[1]
R0.summary$Tmaxh[(R0.summary$strain == "NY10") & (R0.summary$params == "inf")] <- Tmax.CI[2]


R0.summary


# View(R0.summary)
