
using LinearAlgebra
using Turing
using AbstractMCMC
using AbstractGPs, Random, Distributed
using AdvancedHMC
using KernelFunctions, Stheno
using RCall
using Distributions, DistributionsAD, ConditionalDists
using Distances

K = ...;

addprocs(4)

function mat_euclidean(mat)
  N, P = size(mat)
  res = zeros(N, N)
  for i ∈ 1:1:N
    for j ∈ 1:1:(N-1)
      res[i, j] = sum((mat[i, :] - mat[j, :]).^2)
      res[j, i] = res[i, j]
    end
  end
  return res
end

mat_euclidean(zin)


# filldist
xin = Random.rand(100, 4)
yin = Random.rand(100)
zin = Random.rand(100, 3)
K = KernelFunctions.GaussianKernel
a_sigma_in = 2
b_sigma_in = 16
a_lambda_in = 1
b_lambda_in = 0.1
a_pi0_in = 4
b_pi0_in = 2

@model function skeleton_bkmr(x, y, z, K, a_sigma, b_sigma, a_lambda, b_lambda, a_pi0, b_pi0)
  nrow, ncol = size(x)
  
  u = Random.rand(nrow, 1)
  Dz = mat_euclidean(z)
  h = Random.rand(nrow, 1)
  u ~ Uniform(0, 1)
  Σ = u .* K(Dz).metric
  Σ = round.(Σ, sigdigits = 8)
  Σ = Symmetric(Σ)
  h ~ MultivariateNormal(repeat([0.], nrow), Σ)
  # cholesky factorization failed...
  # for i ∈ 1:1:nrow
  y ~ MultivariateNormal(h .+ (x * beta), (tau^-2) .* Diagonal(repeat([1], nrow)))
  # end

  # variable selection
  # Priors
  tau ~ Gamma(a_sigma, b_sigma)
  lambda = u * (tau ^ -2)
  lambda ~ Gamma(a_lambda, b_lambda) # lambda=tau sigma^-2 in the original article; tau is u, sigma^-2 is tau
  pi0 ~ Beta(a_pi0, b_pi0)

  for m ∈ 1:1:ncol
    # P0... density with point mass at 0
    # where pdf is Unif^-1
    r[m] ~ delta[m] .* Distributions.pdf(Uniform(), r[m]) + (1 - delta[m]) * Distributions.pdf(Uniform(), 0)
    delta[m] ~ Bernoulli(pi0)
  end

  for g ∈ 1:1:ngroups
    Delta[g] ~ Multinomial(omega[g], pdf_G) # g ∈ 1, ..., G
    omega[g] ~ Bernoulli(pi0)
  end

end


model = skeleton_bkmr(xin, yin, zin, K, 
    a_sigma_in, b_sigma_in,
    a_lambda_in, b_lambda_in,
    a_pi0_in, b_pi0_in)
sample(model, Turing.NUTS(200, 0.65), 1000)


@everywhere using Turing
@everywhere using KernelFunctions
@everywhere using Random
@everywhere using Distributions

@everywhere begin
xin = Random.rand(100, 4)
yin = Random.rand(100)
zin = Random.rand(100, 3)
K = KernelFunctions.GaussianKernel()
a_sigma_in = 0.0625
b_sigma_in = 0.5
a_lambda_in = 1
b_lambda_in = 0.1
a_pi0_in = 4
b_pi0_in = 2
end

# Define a model on all processes.
@everywhere function skeleton_bkmr(
    x, y, z, K,
    a_sigma, b_sigma,
    a_lambda, b_lambda,
    a_pi0, b_pi0,
    ngroups = 1)
    nrow, ncol = size(x)
    
    u = 0
    h = Random.rand(nrow, 1)
    u ~ Uniform(-1, 1)
    h ~ MultivariateNormal(zeros(nrow), u .* K((z - z')^2))
    # for i ∈ 1:1:nrow
    y ~ MultivariateNormal(h .+ x * beta, (tau^-2) .* Diagonal(repeat([1], nrow)))
    # end
  
    # variable selection
    # Priors
    tau ~ Gamma(a_sigma, b_sigma)
    lambda = u * (tau ^ -2)
    lambda ~ Gamma(a_lambda, b_lambda) # lambda=tau sigma^-2 in the original article; tau is u, sigma^-2 is tau
    pi0 ~ Beta(a_pi0, b_pi0)

    for m ∈ 1:1:ncol
      # P0... density with point mass at 0
      # where pdf is Unif^-1
      r[m] ~ delta[m] .* Distributions.pdf(Uniform(), r[m]) + (1 - delta[m]) * Distributions.pdf(Uniform(), 0)
      delta[m] ~ Bernoulli(pi0)
    end
  
    for g ∈ 1:1:ngroups
      Delta[g] ~ Multinomial(omega[g], pdf_G) # g ∈ 1, ..., G
      omega[g] ~ Bernoulli(pi0)
    end
  
  end

@everywhere model = skeleton_bkmr(xin, yin, zin, K, 
    a_sigma_in, b_sigma_in,
    a_lambda_in, b_lambda_in,
    a_pi0_in, b_pi0_in)
    
# Sample four chains using multiple processes, each with 1000 samples.
  

#####

@everywhere using Turing

# Define a model on all processes.
@everywhere @model function gdemo(x)
    s² ~ InverseGamma(2, 3)
    m ~ Normal(0, sqrt(s²))

    for i in eachindex(x)
        x[i] ~ Normal(m, sqrt(s²))
    end
end

# Declare the model instance everywhere.
@everywhere model = gdemo([1.5, 2.0])

# Sample four chains using multiple processes, each with 1000 samples.
sample(model, NUTS(), MCMCDistributed(), 1000, 4)


## Ref
# n_data = 50
# gplvm_sparse = GPLVM_sparse(dat[1:n_data, :], ndim)

# chain_gplvm_sparse = sample(gplvm_sparse, NUTS(), 500)