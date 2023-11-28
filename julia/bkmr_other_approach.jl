# with helps from ChatGPT by following prompt:
# Can you give me a Julia Turing code for Bayesian Kernel Machine Regression, where
# one takes an exposure matrix Z, a covariate matrix X, and an outcome matrix Y.
# Exposures are supposed to be transformed with a dose-response function, say, h(x), and
# for every pair of exposures are interrelated each other for which is modeled with a kernel,
# which deemed to result in a Gaussian Process. On top of that, there is a grouping procedure
# which is governed by multinomial distribution
# under choosing probability following a binomial distribution for every exposure. 

using Turing, Distributions, LinearAlgebra, AbstractGPs, KernelFunctions, Stheno
using DynamicHMC
using Turing: Variational


sigmoid(z::Real) = one(z) / (one(z) + exp(-z))

# Define your dose-response function h(x)
function h(x::Real)
    transformed_x = sigmoid(x)
    # Implement the dose-response transformation here
    return transformed_x
end

h(3)
h.([3. 4. 3.2 2. 1.5; 1 3 2 4 0])
# Define the kernel function for the Gaussian Process
function K(x, y)
    # elementwise difference
    dxy = x - y
    KernelFunctions.SqExponentialKernel(dxy)
    # Implement your kernel logic here
    # return kernel_value
end

@model function bkmr(Y, Z, X)
    # Number of observations and exposures
    n = size(Y, 1)
    p = size(Z, 2)
    k = size(X, 2)

    # Priors for Gaussian Process parameters
    σ² ~ InverseGamma(2, 3)
    # l ~ filldist(Gamma(2, 3), p)
    
    # Transform exposures
    Z_transformed = h.(Z)

    # Define the Gaussian Process
    gp = GP(zeros(n), (SqExponentialKernel() ∘ ScaleTransform(σ²)))
    # gp = GP(τ .* (Matern32Kernel() ∘ ScaleTransform(σ²)))

    # minit = reshape(mean(Z_transformed, dims = 1), 1, 5)
    # Gaussian Process for the transformed exposures
    # f = gp(minit, 1E-3)
    delta ~ Bernoulli(0.5)
    # slab-and-spike prior
    r ~ delta * Normal(0, 1) + (1-delta) * 0.3
    f = cov(gp(r .* ColVecs(Z_transformed')))
    fx ~ filldist(MvNormal(zeros(n), f), n)

    # Grouping probabilities (assuming a simple case)
    π ~ Dirichlet(ones(p) .* 3)  # Adjust as needed

    # Covariate effects
    β ~ MvNormal(zeros(k), diagm(ones(k)))
    a ~ Normal(0, 1)
    # Model for the outcome
    for i in 1:n
        μ = a + dot(X[i, :], β) + fx[i]
        Y[i] ~ Normal(μ, σ²)
    end
end

# Example data (replace with your actual data)
Z = rand(100, 5)  # 100 observations, 5 exposures
X = rand(100, 3)  # 100 observations, 3 covariates
Y = rand(100, 1)     # 100 outcomes

# Running the model
h.(Z)

model = bkmr(Y, Z, X)
# chain = sample(model, NUTS(0.66), 1000)

using Zygote, ReverseDiff
Turing.setadbackend(:reversediff)
# VI
# input length mismatch?
viest = vi(model, ADVI(4, 1000))
mfest = Variational.meanfield(model)

opt = Variational.DecayedADAGrad(1e-2, 1.1, 0.9)
advi = ADVI(10, 1000)
qvi = vi(model, advi, mfest; optimizer=opt)
