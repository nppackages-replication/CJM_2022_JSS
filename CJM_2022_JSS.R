#--------------------------------------------------------------------------------
# lpdensity: Local Polynomial Density Estimation and Inference
# Matias D. Cattaneo, Michael Jansson, and Xinwei Ma
# Replication Code
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# Part I: Figures and illustrations
#--------------------------------------------------------------------------------

rm(list=ls())

#--------------------------------------------------------------------------------
# Install and load the "lpdensity" and "ggplot2" packages
#--------------------------------------------------------------------------------

# install.packages("ggplot2")
# install.packages("lpdensity")

library("ggplot2")
library("lpdensity")

#----------------------------------------
# Generate data
#----------------------------------------

set.seed(42)
data <- rnorm(4000, mean = -1)
data <- data[data < 0]
data <- -1 * data[1:2000]

#----------------------------------------
# Figure 1, panel (a)
#----------------------------------------

dataHist <- as.data.frame(data)
colnames(dataHist) <- c("v1")
dataHist$pdf <- dnorm(dataHist$v1, mean = 1, sd = 1) / pnorm(0, mean = 1, sd = 1, lower.tail = FALSE)
ggplot() + geom_histogram(data=dataHist, aes(x=v1, y=..density..), breaks=seq(0, 4, 0.2), fill=2, col="white", alpha=0.6) +
  theme_bw() + labs(x = "", y = "Density") +
  geom_line(data=dataHist, aes(x=v1, y=pdf))

#----------------------------------------
# Figure 1, panel (b)
#----------------------------------------

model2 <- lpdensity(data, bw = 0.5, grid = seq(0, 4, 0.05))
plot(model2) + theme(legend.position = "none")

#----------------------------------------
# lpdensity(): Estimation with bandwidth 0.5 on provided grid points
#----------------------------------------

model1 <- lpdensity(data, bw = 0.5, grid = seq(0, 4, 0.5))
summary(model1)

#----------------------------------------
# lpdensity(): extracting estimation
#   results
#----------------------------------------

model1$Estimate

#----------------------------------------
# lpdensity(): conventional inference
#----------------------------------------

summary(lpdensity(data, bw = 0.5, p = 2, q = 2))

#----------------------------------------
# lpdensity(): customizing screen output
#----------------------------------------

set.seed(123) # fix the random seed for critical value simulation
summary(model1, alpha = 0.01, sep = 3, grid = c(0, 0.5, 1, 2), CIuniform = TRUE)

#----------------------------------------
# lpdensity(): inconsistent density
#   estimation using partial sample
#----------------------------------------

lpdensity(data[data < 1.5], bw = 0.5, grid = 1.5)$Estimate[, "f_p"]
lpdensity(data[data > 1.5], bw = 0.5, grid = 1.5)$Estimate[, "f_p"]
dnorm(1.5, mean = 1, sd = 1) / pnorm(0, mean = 1, sd = 1, lower.tail = FALSE) # true density at 1.5

#----------------------------------------
# lpdensity(): consistent density
#   estimation using partial sample and
#   option "scale"
#----------------------------------------

lpdensity(data[data < 1.5], bw = 0.5, grid = 1.5,
          scale = sum(data < 1.5)/2000)$Estimate[, "f_p"]
lpdensity(data[data > 1.5], bw = 0.5, grid = 1.5,
          scale = sum(data > 1.5)/2000)$Estimate[, "f_p"]

#----------------------------------------
# plot(): customization
#----------------------------------------

plot(model2, CItype="line") + theme(legend.position = "none")
plot(model2, type="points", CItype="ebar", grid = seq(0, 4, 0.5)) +
  theme(legend.position = "none")
plot(model2, hist = TRUE, histData = data, histBreaks = seq(0, 4, 0.2)) +
  theme(legend.position = "none")
set.seed(123) # fix the random seed for critical value simulation
plot(model2, alpha=0.1, CIuniform = TRUE) +
  theme(legend.position = "none")

#----------------------------------------
# lpbwdensity(): illustration
#----------------------------------------

model1bw <- lpbwdensity(data, grid = seq(0, 4, 0.5))
summary(model1bw)

#----------------------------------------
# lpdensity(): automatic bandwidth
#   selection
#----------------------------------------

model5 <- lpdensity(data, grid = seq(0, 4, 0.5), bwselect = "imse-dpi")
summary(model5)

#----------------------------------------
# lpdensity(): undersmoothing
#----------------------------------------

# Estimation and plot using IMSE bandwidth
model6bwIMSE <- lpbwdensity(data, grid = seq(0, 4, 0.05),
                            bwselect = "imse-dpi")
model6 <- lpdensity(data, grid = seq(0, 4, 0.05),
                    bw = model6bwIMSE$BW[, "bw"])
plot(model6) + theme(legend.position = "none")

# Estimation and plot using half IMSE bandwidth
model7 <- lpdensity(data, grid = seq(0, 4, 0.05),
                    bw = model6bwIMSE$BW[, "bw"] / 2)
plot(model7) + theme(legend.position = "none")

#--------------------------------------------------------------------------------
# Part II: Simulation - Truncated Normal
#--------------------------------------------------------------------------------

rm(list=ls())

#--------------------------------------------------------------------------------
# Install and load packages
#--------------------------------------------------------------------------------

# install.packages("lpdensity")
# install.packages("KernSmooth")
# install.packages("ks")
# install.packages("np")
# install.packages("nprobust")
# install.packages("plugdensity")

library("lpdensity")
library("KernSmooth")
library("ks")
library("np")
library("nprobust")
library("plugdensity")

#----------------------------------------
# Simulation setup
#----------------------------------------

repe <- 2000 # number of simulations

n <- 1000 # sample size
xList <- c(-1.5, -0.2, 0) # evaluation points

#----------------------------------------
# Data generating process
#----------------------------------------

# N(-1, 1) truncated from above at 0
dataGen <- function(n) {
  temp <- c()
  while (length(temp) < n) {
    temp <- c(temp, rnorm(n, mean = -1, sd = 1))
    temp <- temp[temp <= 0]
  }
  return(temp[1:n])
}

#----------------------------------------
# Simulation
#----------------------------------------
ptm <- proc.time()

for (args in 1:length(xList)) {
  # Seed
  set.seed(42)
  x <- xList[args]

  cat(paste("The job ID is",   args, "\n", sep=" "))
  cat(paste("x is",               x, "\n", sep=" "))

  resultMatrix <- matrix(NA, nrow = repe, ncol = 44) # collecting results
  truePar <- dnorm(x, mean=-1, sd=1) / pnorm(0, mean=-1, sd=1) # true parameter

  interior <- (args == 1) # interior

  for (i in 1:repe) {
    cat(i); cat("\n")
    data <- dataGen(n)

    # KernSmooth::bkde
    if (interior) {
      temp <- KernSmooth::bkde(x = data)
      resultMatrix[i, 1] <- NA
      tempindex <- which.min(abs(temp$x - x))
      resultMatrix[i, 2] <- temp$y[tempindex]  - truePar
      resultMatrix[i, 3] <- NA
      resultMatrix[i, 4] <- NA
    }

    # KernSmooth::locpoly
    resultMatrix[i, 5] <- KernSmooth::dpik(x = data)
    temp <- KernSmooth::locpoly(x = data, drv = 0, degree = 1, bandwidth = resultMatrix[i, 5])
    tempindex <- which.min(abs(temp$x - x))
    resultMatrix[i, 6] <- temp$y[tempindex]  - truePar
    resultMatrix[i, 7] <- NA
    resultMatrix[i, 8] <- NA

    # ks::kdde
    if (interior) {
      temp <- ks::kdde(x = data, deriv.order = 0, eval.points = x)
      resultMatrix[i, 9 ] <- temp$h
      resultMatrix[i, 10] <- temp$estimate  - truePar
      resultMatrix[i, 11] <- NA
      resultMatrix[i, 12] <- NA
    }

    # ks::kde
    if (interior) {
      temp <- ks::kde(x = data, eval.points = x)
      resultMatrix[i, 13] <- temp$h
      resultMatrix[i, 14] <- temp$estimate  - truePar
      resultMatrix[i, 15] <- NA
      resultMatrix[i, 16] <- NA
    }

    # np::npudens
    if (interior) {
      temp <- np::npudens(tdat = data, edat = x)
      resultMatrix[i, 17] <- temp$bw
      resultMatrix[i, 18] <- temp$dens  - truePar
      resultMatrix[i, 19] <- 1 * (abs(temp$dens  - truePar) / temp$derr <= 1.96)
      resultMatrix[i, 20] <- 2 * 1.96 * temp$derr
    }

    # np::npuniden.boundary
    temp <- np::npuniden.boundary(X = data, a = min(data), b = max(data))
    tempindex <- which.min(abs(data - x))
    resultMatrix[i, 21] <- temp$h
    resultMatrix[i, 22] <- temp$f[tempindex] - truePar
    resultMatrix[i, 23] <- 1 * (abs(temp$f[tempindex] - truePar) / temp$sd.f[tempindex] <= 1.96)
    resultMatrix[i, 24] <- 2 * 1.96 * temp$sd.f[tempindex]

    # nprobust::kdrobust
    if (interior) {
      temp <- nprobust::kdrobust(x = data, bwselect = "imse-dpi")
      temp <- nprobust::kdrobust(x = data, eval = x, h = temp$Estimate[1, 2])
      resultMatrix[i, 25] <- temp$Estimate[1, 2]
      resultMatrix[i, 26] <- temp$Estimate[1, 5] - truePar
      resultMatrix[i, 27] <- 1 * (abs(temp$Estimate[1, 6] - truePar) / temp$Estimate[1, 8] <= 1.96)
      resultMatrix[i, 28] <- 2 * 1.96 * temp$Estimate[1, 8]
    }

    # plugdensity::plugin.density
    if (interior) {
      temp <- plugdensity::plugin.density(x = data, xout = x)
      resultMatrix[i, 29] <- temp$bw
      resultMatrix[i, 30] <- temp$y  - truePar
      resultMatrix[i, 31] <- NA
      resultMatrix[i, 32] <- NA
    }

    # stats::density
    if (interior) {
      temp <- stats::density(x = data, n = 1, from = x, to = x)
      resultMatrix[i, 33] <- temp$bw
      resultMatrix[i, 34] <- temp$y - truePar
      resultMatrix[i, 35] <- NA
      resultMatrix[i, 36] <- NA
    }

    # lpdensity::lpdensity, MSE bw
    temp <- lpdensity(data = data, grid = x, bwselect = "mse-dpi", p = 2, v = 1, kernel = "triangular")
    resultMatrix[i, 37] <- temp$Estimate[1, "bw"]
    resultMatrix[i, 38] <- temp$Estimate[1, "f_p"] - truePar
    resultMatrix[i, 39] <- 1 * (abs(temp$Estimate[1, "f_q"] - truePar) / temp$Estimate[1, "se_q"] <= 1.96)
    resultMatrix[i, 40] <- 2 * 1.96 * temp$Estimate[1, "se_q"]

    # lpdensity::lpdensity, IMSE bw
    temp <- lpdensity(data = data, bwselect = "imse-dpi", p = 2, v = 1, kernel = "triangular")
    temp <- lpdensity(data = data, grid = x, bw = temp$Estimate[1, 2], p = 2, v = 1, kernel = "triangular")
    resultMatrix[i, 41] <- temp$Estimate[1, "bw"]
    resultMatrix[i, 42] <- temp$Estimate[1, "f_p"] - truePar
    resultMatrix[i, 43] <- 1 * (abs(temp$Estimate[1, "f_q"] - truePar) / temp$Estimate[1, "se_q"] <= 1.96)
    resultMatrix[i, 44] <- 2 * 1.96 * temp$Estimate[1, "se_q"]
  }

  fname <- paste("ResultNormal", "x", args, ".txt", sep = "")
  write.table(resultMatrix, file = fname, sep = ",", row.names = FALSE, col.names = FALSE)
}

proc.time() - ptm

#--------------------------------------------------------------------------------
# Part II: Simulation - Exponential
#--------------------------------------------------------------------------------

#----------------------------------------
# Data generating process
#----------------------------------------

# Exp(1) distribution
dataGen <- function(n) {
  return(rexp(n, rate = 1))
}

#----------------------------------------
# Simulation
#----------------------------------------
ptm <- proc.time()

for (args in 1:length(xList)) {
  # Seed
  set.seed(42)
  x <- xList[args]

  cat(paste("The job ID is",   args, "\n", sep=" "))
  cat(paste("x is",               x, "\n", sep=" "))

  resultMatrix <- matrix(NA, nrow = repe, ncol = 44) # collecting results
  truePar <- dexp(x, rate = 1) # true parameter

  interior <- (args == 1) # interior

  for (i in 1:repe) {
    cat(i); cat("\n")
    data <- dataGen(n)

    # KernSmooth::bkde
    if (interior) {
      temp <- KernSmooth::bkde(x = data)
      resultMatrix[i, 1] <- NA
      tempindex <- which.min(abs(temp$x - x))
      resultMatrix[i, 2] <- temp$y[tempindex]  - truePar
      resultMatrix[i, 3] <- NA
      resultMatrix[i, 4] <- NA
    }


    # KernSmooth::locpoly
    resultMatrix[i, 5] <- KernSmooth::dpik(x = data)
    temp <- KernSmooth::locpoly(x = data, drv = 0, degree = 1, bandwidth = resultMatrix[i, 5])
    tempindex <- which.min(abs(temp$x - x))
    resultMatrix[i, 6] <- temp$y[tempindex]  - truePar
    resultMatrix[i, 7] <- NA
    resultMatrix[i, 8] <- NA

    # ks::kdde
    if (interior) {
      temp <- ks::kdde(x = data, deriv.order = 0, eval.points = x)
      resultMatrix[i, 9 ] <- temp$h
      resultMatrix[i, 10] <- temp$estimate  - truePar
      resultMatrix[i, 11] <- NA
      resultMatrix[i, 12] <- NA
    }

    # ks::kde
    if (interior) {
      temp <- ks::kde(x = data, eval.points = x)
      resultMatrix[i, 13] <- temp$h
      resultMatrix[i, 14] <- temp$estimate  - truePar
      resultMatrix[i, 15] <- NA
      resultMatrix[i, 16] <- NA
    }

    # np::npudens
    if (interior) {
      temp <- np::npudens(tdat = data, edat = x)
      resultMatrix[i, 17] <- temp$bw
      resultMatrix[i, 18] <- temp$dens  - truePar
      resultMatrix[i, 19] <- 1 * (abs(temp$dens  - truePar) / temp$derr <= 1.96)
      resultMatrix[i, 20] <- 2 * 1.96 * temp$derr
    }

    # np::npuniden.boundary
    temp <- np::npuniden.boundary(X = data, a = min(data), b = max(data))
    tempindex <- which.min(abs(data - x))
    resultMatrix[i, 21] <- temp$h
    resultMatrix[i, 22] <- temp$f[tempindex] - truePar
    resultMatrix[i, 23] <- 1 * (abs(temp$f[tempindex] - truePar) / temp$sd.f[tempindex] <= 1.96)
    resultMatrix[i, 24] <- 2 * 1.96 * temp$sd.f[tempindex]

    # nprobust::kdrobust
    if (interior) {
      temp <- nprobust::kdrobust(x = data, bwselect = "imse-dpi")
      temp <- nprobust::kdrobust(x = data, eval = x, h = temp$Estimate[1, 2])
      resultMatrix[i, 25] <- temp$Estimate[1, 2]
      resultMatrix[i, 26] <- temp$Estimate[1, 5] - truePar
      resultMatrix[i, 27] <- 1 * (abs(temp$Estimate[1, 6] - truePar) / temp$Estimate[1, 8] <= 1.96)
      resultMatrix[i, 28] <- 2 * 1.96 * temp$Estimate[1, 8]
    }

    # plugdensity::plugin.density
    if (interior) {
      temp <- plugdensity::plugin.density(x = data, xout = x)
      resultMatrix[i, 29] <- temp$bw
      resultMatrix[i, 30] <- temp$y  - truePar
      resultMatrix[i, 31] <- NA
      resultMatrix[i, 32] <- NA
    }

    # stats::density
    if (interior) {
      temp <- stats::density(x = data, n = 1, from = x, to = x)
      resultMatrix[i, 33] <- temp$bw
      resultMatrix[i, 34] <- temp$y - truePar
      resultMatrix[i, 35] <- NA
      resultMatrix[i, 36] <- NA
    }

    # lpdensity::lpdensity, MSE bw
    temp <- lpdensity(data = data, grid = x, bwselect = "mse-dpi", p = 2, v = 1, kernel = "triangular")
    resultMatrix[i, 37] <- temp$Estimate[1, "bw"]
    resultMatrix[i, 38] <- temp$Estimate[1, "f_p"] - truePar
    resultMatrix[i, 39] <- 1 * (abs(temp$Estimate[1, "f_q"] - truePar) / temp$Estimate[1, "se_q"] <= 1.96)
    resultMatrix[i, 40] <- 2 * 1.96 * temp$Estimate[1, "se_q"]

    # lpdensity::lpdensity, IMSE bw
    temp <- lpdensity(data = data, bwselect = "imse-dpi", p = 2, v = 1, kernel = "triangular")
    temp <- lpdensity(data = data, grid = x, bw = temp$Estimate[1, 2], p = 2, v = 1, kernel = "triangular")
    resultMatrix[i, 41] <- temp$Estimate[1, "bw"]
    resultMatrix[i, 42] <- temp$Estimate[1, "f_p"] - truePar
    resultMatrix[i, 43] <- 1 * (abs(temp$Estimate[1, "f_q"] - truePar) / temp$Estimate[1, "se_q"] <= 1.96)
    resultMatrix[i, 44] <- 2 * 1.96 * temp$Estimate[1, "se_q"]
  }

  fname <- paste("ResultExp", "x", args, ".txt", sep = "")
  write.table(resultMatrix, file = fname, sep = ",", row.names = FALSE, col.names = FALSE)
}

proc.time() - ptm

#--------------------------------------------------------------------------------
# Part II: Simulation - Generate Table 2
#--------------------------------------------------------------------------------

xList <- c(1.5, 0.2, 0) # evaluation point
nModel <- 11 # number of models

Output <- matrix(NA, ncol=12, nrow=nModel*3)
for (x in 1:length(xList)) {
  fname <- paste("ResultNormal", "x", x, ".txt", sep = "")
  temp <- read.table(fname, sep=",")

  for (j in 1:nModel) {
    index <- (j-1) * 4
    Output[(x-1)*11 + j, 1] <- mean(temp[, c(index+1)])
    Output[(x-1)*11 + j, 2] <- abs(mean(temp[, c(index+2)]))
    Output[(x-1)*11 + j, 3] <- sd(temp[, c(index+2)])
    Output[(x-1)*11 + j, 4] <- sqrt(abs(mean(temp[, c(index+2)]))^2 + sd(temp[, c(index+2)])^2)
    Output[(x-1)*11 + j, 5] <- mean(temp[, c(index+3)])
    Output[(x-1)*11 + j, 6] <- mean(temp[, c(index+4)])
  }

  fname <- paste("ResultExp", "x", x, ".txt", sep = "")
  temp <- read.table(fname, sep=",")

  for (j in 1:nModel) {
    index <- (j-1) * 4
    Output[(x-1)*11 + j, 6 + 1] <- mean(temp[, c(index+1)])
    Output[(x-1)*11 + j, 6 + 2] <- abs(mean(temp[, c(index+2)]))
    Output[(x-1)*11 + j, 6 + 3] <- sd(temp[, c(index+2)])
    Output[(x-1)*11 + j, 6 + 4] <- sqrt(abs(mean(temp[, c(index+2)]))^2 + sd(temp[, c(index+2)])^2)
    Output[(x-1)*11 + j, 6 + 5] <- mean(temp[, c(index+3)])
    Output[(x-1)*11 + j, 6 + 6] <- mean(temp[, c(index+4)])
  }
}

Output <- Output[c(1:11, 13, 17, 21, 22, 24, 28, 32, 33), ]
round(Output, 3)
