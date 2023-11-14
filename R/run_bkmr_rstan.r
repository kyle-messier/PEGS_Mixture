library(rstan)

options(mc.cores = 8)
rstan_options(auto_write = TRUE, threads_per_chain = 1)

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
zin <- matrix(runif(100 * 3), nrow = 100)
K <- "GaussianKernel"  # Define the Gaussian Kernel function
a_sigma_in <- 2
b_sigma_in <- 16
a_lambda_in <- 1
b_lambda_in <- 0.1
a_pi0_in <- 4
b_pi0_in <- 2

# Define the Stan model
stan_model <- "
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
