using LinearAlgebra
using Turing
using AbstractGPs, Random, Distributed
using AdvancedHMC
using KernelFunctions, Stheno
using RCall
using Distributions, DistributionsAD, ConditionalDists

K = ...;

addprocs(4)

@model function skeleton_bkmr(x, y, z, K, a_sigma, b_sigma, a_lambda, b_lambda, a_pi0, b_pi0)
  nrow, ncol = size(x)
  
  # Priors
  tau ~ Gamma(a_sigma, b_sigma)
  lambda ~ Gamma(a_lambda, b_lambda) # lambda=tau sigma^-2 in the original article; tau is u, sigma^-2 is tau
  pi0 ~ Beta(a_pi0, b_pi0)

  # for i ∈ 1:1:nrow
  y ~ Normal(h .+ x * beta, (tau^-2) .* Diagonal(repeat([1], nrow)))
  # end
  h ~ MultivariateNormal(zeros(nrow), u .* z .* K)

  # variable selection
  
  for m ∈ ncomps
    # P0... density with point mass at 0
    r[m] ~ delta[m] .* pdf(r[m]) + (1 - delta[m]) * P0
    delta[m] ~ Bernoulli(pi0)
  end

  for g ∈ ngroups
    Delta ~ Multinomial(omega[g], pdf_G) # g ∈ 1, ..., G
    omega[g] ~ Bernoulli(pi0)
  end

end

# filldist
xin = Random.rand(100, 4)
yin = Random.rand(100)
zin = Random.rand(100, 3)
K = KernelFunctions.GaussianKernel()
a_sigma_in = 2
b_sigma_in = 16
a_lambda_in = 1
b_lambda_in = 0.1
a_pi0_in = 4
b_pi0_in = 2

@everywhere using Turing

# Define a model on all processes.
@everywhere function skeleton_bkmr(
    x, y, z, K,
    a_sigma, b_sigma,
    a_lambda, b_lambda,
    a_pi0, b_pi0)
    nrow, ncol = size(x)
    
    # Priors
    tau ~ Gamma(a_sigma, b_sigma)
    lambda ~ Gamma(a_lambda, b_lambda) # lambda=tau sigma^-2 in the original article; tau is u, sigma^-2 is tau
    pi0 ~ Beta(a_pi0, b_pi0)
  
    # for i ∈ 1:1:nrow
    y ~ Normal(h .+ x * beta, (tau^-2) .* Diagonal(repeat([1], nrow)))
    # end
    h ~ MultivariateNormal(zeros(nrow), u .* K((k - k')^2))
  
    # variable selection
    
    for m ∈ ncomps
      # P0... density with point mass at 0
      # where pdf is Unif^-1
      r[m] ~ delta[m] .* pdf(r[m]) + (1 - delta[m]) * P0
      delta[m] ~ Bernoulli(pi0)
    end
  
    for g ∈ ngroups
      Delta[g] ~ Multinomial(omega[g], pdf_G) # g ∈ 1, ..., G
      omega[g] ~ Bernoulli(pi0)
    end
  
  end

@everywhere model = skeleton_bkmr(xin, yin, zin, K, 
    a_sigma_in, b_sigma_in,
    a_lambda_in, b_lambda_in,
    a_pi0_in, b_pi0_in)
    
# Sample four chains using multiple processes, each with 1000 samples.
sample(model, NUTS(), MCMCDistributed(), 1000, 4)
  

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