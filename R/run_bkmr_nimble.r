if (dir.exists("~/r-libs")) {
  .libPaths("~/r-libs")
}

library(nimble)
# nimbleOptions('showCompilerOutput' = TRUE)
Sys.setenv(LD_LIBRARY_PATH = paste("/ddn/gs1/home/songi2/nimble-dylib", Sys.getenv("LD_LIBRARY_PATH"), sep = ":"))
Sys.getenv("LD_LIBRARY_PATH")

# Function to calculate the Euclidean distance matrix
mat_euclidean <- function(mat) {
  N <- nrow(mat)
  res <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:(N - 1)) {
      res[i, j] <- sum((mat[i, ] - mat[j, ])^2)
      res[j, i] <- res[i, j]
    }
  }
  return(res)
}

# Generate sample data
set.seed(123)
xin <- matrix(runif(100 * 4), nrow = 100)
yin <- runif(100)
zin <- matrix(runif(100 * 8), nrow = 100)
N <- nrow(xin)
K <- "GaussianKernel"  # Define the Gaussian Kernel function
a_sigma_in <- 2
b_sigma_in <- 16
a_lambda_in <- 1
b_lambda_in <- 0.1
a_pi0_in <- 4
b_pi0_in <- 2
Pi_s_in <- c(0.25, 0.3, 0.33, 0.1)

data_in <-
  list(x = xin, y = yin, z = zin)
init_in <-
  list(a_pi0 = a_pi0_in, b_pi0 = b_pi0_in, tau = 3, Pi_s = Pi_s_in,
  Sigma = cor(zin), tau=0.01, beta = c(1,1,1,1), prob_omega = 0.5,
  mu = rep(0, nrow(xin)))
const_in <- list(N = nrow(xin), P = ncol(xin), ngroups = 2L,
  h = rep(0, 100))

mat_euclidean <- nimbleFunction(
  run = function(mat = double(2)) {
  N <- dim(mat)[1]
  res <- matrix(0, N, N)
  for (i in 1:N) {
    for (j in 1:(N - 1)) {
      res[i, j] <- sum((mat[i, ] - mat[j, ])^2)
      res[j, i] <- res[i, j]
    }
  }
  return(res)
  returnType(double(2))
}
)
Cmat_euclidean <- compileNimble(mat_euclidean)

# Define the nimble model
nimble_model <- nimbleCode({
  #inith <- rep(0, N)
  #hprec <- Sigma

  for (om in 1:ngroups) {
    omega[om] ~ dbern(prob_omega)
  }
  for (g in 1:ngroups) {
    Delta[g] ~ ddirch(omega[g])
  }  
  
  u ~ dunif(0, 1)
  #e ~ dnorm(0, 1)
  Sigma <- u * mat_euclidean(x)
  for (j in 1:P) {
    beta[j] ~ dunif(-1, 1)
  }
  tau ~ dgamma(1, 0.0005)
  a_pi0 ~ dunif(0, 1)
  b_pi0 ~ dunif(0, 1)
  pi0 ~ dbeta(a_pi0, b_pi0)
  for (d in 1:P) {
    delta[d] ~ dbern(pi0)
  }
  
  for (k in 1:N) {
    #h[k] ~ dnorm(0, 1)
    #for (p in 1:P) {
    mu[k] <- hx[k] + inprod(beta[1:P], x[k,])
    #}
  }
  yprec[1:N,1:N] <- (tau^-2) * diag(rep(1, N))
  y[1:N] ~ dmnorm(mu[1:N], Sigma[1:N, 1:N])

})

# nimble_model_code <- nimble::nimbleCode(nimble_model)

nimble_fit <- nimble::nimbleModel(code = nimble_model,
  data = data_in,
  inits = init_in,
  constants = const_in,
  check = TRUE)

nimble_mcmcrun <- nimble::nimbleMCMC(
  model = nimble_fit,
  niter = 1e4L,
  nchains = 4,
  nburnin = 3e3L,
  WAIC = TRUE,
  summary = TRUE)
# nimble_model <- nimbleCode({
#   u ~ dunif(0,1)
#   h ~ dmnorm(rep(0, nrow), Sigma)
#   y ~ dmnorm(h + x * beta, (tau^-2) * diag(rep(1, nrow)))

#   pi0 ~ dbeta(a_pi0, b_pi0)
#   delta ~ dbern(pi0)
#   for (om in 1:ngroups) {
#     omega[om] ~ dbern(0.5)
#   }
#   for (g in 1:ngroups) {
#     Delta[g] ~ ddirch(omega[g])
#   }
# })


nimble_model <- "
functions {
  matrix mat_euclidean(matrix mat) {
    int N = rows(mat);
    matrix[N, N] res;
    for (i in 1:N) {
      for (j in 1:(N - 1)) {
        res[i, j] = sum((mat[i, ] - mat[j, ])^2);
        res[j, i] = res[i, j];
      }
    }
    return res;
  }
}

data {
  int nrow;
  int ncol;
  int ncolz;
  int ngroups;
  matrix[nrow, ncol] x;
  vector[nrow] y;
  matrix[nrow, ncolz] z;
  // real pi0;
  real a_sigma;
  real b_sigma;
  real a_lambda;
  real b_lambda;
  real a_pi0;
  real b_pi0;
}

transformed data {
  matrix[nrow, nrow] Dz;
  matrix[nrow, nrow] Sigma;
  real<lower = 0> tau;

  tau = 1 / sqrt(a_sigma / b_sigma);
  Dz = mat_euclidean(z);
  Sigma = Dz * (a_lambda / b_lambda);
}

parameters {
  vector[nrow] u;
  vector[nrow] h;
  vector[ncol] beta;
  // real<lower = 0> tau;
  // integer<lower=1> ngroups;
  vector<lower = 0, upper = 1>[ncol] delta;
  vector<lower = 0, upper = 1>[ngroups] omega;
  vector<lower = 0>[ngroups] Delta;
  real<lower=0, upper=1> pi0;
}

model {
  u ~ uniform(0, 1);
  h ~ multi_normal(rep_vector(0, nrow), Sigma);
  y ~ multi_normal(h + x * beta, (tau^-2) * diag_matrix(rep_vector(1, nrow)));

  //for (m in 1:ncol) {
  // pi0 ~ beta(a_pi0, b_pi0);
  // delta ~ bernoulli(pi0);
  //}
  delta ~ bernoulli(0.5);
  for (om in 1:ngroups)
    omega[om] ~ bernoulli(0.5);

  for (g in 1:ngroups)
    Delta[g] ~ dirichlet(omega[g]);

}

generated quantities {
  // Additional quantities if needed
}

"


# Define the data
data_list <- list(
  nrow = nrow(xin),
  ncol = ncol(xin),
  ncolz = ncol(zin),
  pi0 = 0.5,
  ngroups = 2,
  x = xin,
  y = yin,
  z = zin,
  a_sigma = a_sigma_in,
  b_sigma = b_sigma_in,
  a_lambda = a_lambda_in,
  b_lambda = b_lambda_in,
  a_pi0 = a_pi0_in,
  b_pi0 = b_pi0_in
)

# Compile the model
stan_model <- stan_model(model_code = stan_model)

# Sample from the model
stan_samples <- sampling(stan_model, data = data_list, chains = 4, iter = 1000, warmup = 200, thin = 1)


# example
zinf <- "
data {
  int<lower=0> N;
  int<lower=0> y[N];
}
parameters {
  real<lower=0, upper=1> theta;
  real<lower=0> lambda;
}
model {
  for (n in 1:N) {
    if (y[n] == 0)
      target += log_sum_exp(bernoulli_lpmf(1 | theta),
                            bernoulli_lpmf(0 | theta)
                              + poisson_lpmf(y[n] | lambda));
    else
      target += bernoulli_lpmf(0 | theta)
                  + poisson_lpmf(y[n] | lambda);
  }
}"

nn <- list(N = 1000, y = rnbinom(1000, 0.5, 0.1))

# Compile the model
stan_model <- stan_model(model_code = zinf)

# Sample from the model
stan_samples <- sampling(stan_model, data = nn, chains = 4, iter = 10000, warmup = 3000, thin = 1)

bayesplot::mcmc_hist_by_chain(stan_samples)




####

library(NLinteraction)
set.seed(2024)
n = 1000
p = 5
pc = 3

X = matrix(rnorm(n*p), n, p)

C = matrix(rnorm(n*pc), nrow=n)
bC = matrix(runif(pc, -3, 3), ncol = 1)

TrueH = function(X) {
  return(1.5*(exp(X[,2])*log(4+X[,3])) - 1.6*(X[,4]^2 * X[,5]))
}

Y = 5 + crossprod(t(C), bC) + TrueH(X) + rnorm(n)

NLmod2 = NLint(Y=Y, X=X, C=C, nIter=2000, nBurn=1000, thin=2, nChains=4, ns=4)
NLmod3 = NLint(Y=Y, X=X, C=C, nIter=2000, nBurn=1000, thin=2, nChains=4, ns=2)

NLinteraction::plotSurface2d(NLmod2, X, C, 2, 3)
