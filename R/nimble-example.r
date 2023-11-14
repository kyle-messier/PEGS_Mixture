dat <- read.csv('https://r-nimble.org/nimbleExamples/avandia.csv')
head(dat)

library(nimble, warn.conflicts = FALSE)

x <- dat$controlMI
n <- dat$nControl
y <- dat$avandiaMI
m <- dat$nAvandia

nStudies <- nrow(dat)
data <- list(x = x, y = y)
constants = list(n = n, m = m, nStudies = nStudies)

codeParam <- nimbleCode({
    for(i in 1:nStudies) {
        y[i] ~ dbin(size = m[i], prob = q[i]) # avandia MIs
        x[i] ~ dbin(size = n[i], prob = p[i]) # control MIs
        q[i] <- expit(theta + gamma[i])       # Avandia log-odds
        p[i] <- expit(gamma[i])               # control log-odds
        gamma[i] ~ dnorm(mu, var = tau2)      # study effects
    }
    theta ~ dflat()        # effect of Avandia
    # random effects hyperparameters
    mu ~ dnorm(0, 10)
    tau2 ~ dinvgamma(2, 1)
})

set.seed(9)
inits = list(theta = 0, mu = 0, tau2 = 1, gamma = rnorm(nStudies))

samples <- nimbleMCMC(code = codeParam, data = data, inits = inits,
                      constants = constants, monitors = c("mu", "tau2", "theta", "gamma"),
                      thin = 10, niter = 22000, nburnin = 2000, nchains = 1, 
                      setSeed = TRUE)


par(mfrow = c(1, 4), cex = 1.1, mgp = c(1.8,.7,0))
ts.plot(samples[ , 'theta'], xlab = 'iteration', ylab = expression(theta),
    main = expression(paste('traceplot for ', theta)))
hist(samples[ , 'theta'], xlab = expression(theta), main = 'effect of Avandia')

gammaCols <- grep('gamma', colnames(samples))
gammaMn <- colMeans(samples[ , gammaCols])
hist(gammaMn, xlab = 'posterior means', main = 'random effects distribution')
hist(samples[1000, gammaCols], xlab = 'single draw',
                   main = 'random effects distribution')