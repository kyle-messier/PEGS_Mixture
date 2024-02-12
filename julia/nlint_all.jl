function ActivePred(zeta, k)
    active = Dict{Int,Vector{String}}()
    
    for h in 1:k
        wv = findall(x->x==1, zeta[:,h])
        active[h] = [""]
        if length(wv) > 0
            push!(active[h], join(wv))
        end
        if length(wv) > 1
            for l in 2:length(wv)
                comb = collect(combinations(wv, l)) 
                for j in eachindex(comb)
                    push!(active[h], join(sort(comb[j]))) 
                end
            end
        end
    end
    
    for h in 2:k
        wv = findall(x->x==1, zeta[:,h]) 
        active[h] = [0]
        if length(wv) > 0
            for l in 1:length(wv)
                if sum(zeta[wv[l], 1:(h-1)]) == 0 
                    push!(active[h], wv[l])
                end 
            end
        end
        if length(wv) > 1
            for l in 2:length(wv)
                comb = collect(combinations(wv, l))
                for j in eachindex(comb)
                    strcomb = join(sort(comb[j]))
                    if !(strcomb in vcat(values(active)[1:(h-1)]...))
                        push!(active[h], strcomb) 
                    end
                end 
            end 
        end
    end
    
    return active
end


function updateBetaOne(Y, tempZeta, f_jhi_nc, betaC, sigmaP, tau, k, sigB, 
    Xstar, tempBeta, h, designC, ns, groups, intMax)

    n = size(designC,1)

    activeZ = activePred(tempZeta, k)
    numZero = sum(sum(tempZeta[:,2:end], dims=2) .== 0)

    tempY = Y - sum(f_jhi_nc[:,2:end], dims=2) - designC * betaC

    tempZeta10 = copy(tempZeta) 
    tempZeta00 = copy(tempZeta)
    tempZeta01 = copy(tempZeta)
    tempZeta10[groups,h] .= 1
    tempZeta00[groups,h] .= 0

    wv = findall(x -> x==1, tempZeta00[:,h])
    sizewv = length(wv)

    p00 = p10 = -Inf 
    p01 = fill(-Inf, 2^length(wv) - 1)

    # Calculate probability for reduced model
    if sizewv == 0
        p00 = log(1 - tau[h])
    elseif sizewv == 1 && sum(tempZeta[wv,1:h-1]) == 0
        Xtemp = Xstar[:,wv,:]
        Σ = Diagonal(sigmaP .* sigB .* ones(size(Xtemp,2)))
        μ = zeros(size(Xtemp,2))
        μβ = (Xtemp'Xtemp/sigmaP + inv(Σ)) \ (Xtemp'Y/sigmaP + inv(Σ)*μ) 
        Σβ = inv(Xtemp'Xtemp/sigmaP + inv(Σ))
        p00 = log(tau[h] * (1 - tau[h])) +
            mvn_logpdf(zeros(ns), μ[2:end], Σ[2:end,2:end]) -
            mvn_logpdf(zeros(ns), μβ[2:end], Σβ[2:end,2:end])
    end

    # Calculate probability for full interactions
    if (numZero > 0 && sum(tempZeta10[:,h]) > 1 && sizewv < intMax) || 
        sum(tempZeta10[:,h]) == 1
    # code to calculate p10
    p10 = #...
    
    end

    PossMatTemp = vec(collect(Iterators.product([0,1], length(wv))))
    PossMat = PossMatTemp[1:end-1,:]

    wh01 = findfirst(x-> sum(tempZeta[:,x]) == 0 && x != h)

    # Now calculate probability for reduced model + effect elsewhere 
    if numZero > 0 && #...
    # code to calculate p01
    end

    # Sample 
    maxlog = maximum(p00, p10, p01)  
    updated_p = exp.((p00, p10, p01) .- maxlog) ./ sum(exp.((p00, p10, p01) .- maxlog))
    samp = rand(Categorical(updated_p))

    # Update tempZeta based on samp
    #...

    # Update tempBeta and f_jhi_nc
    #...

    return (f_jhi_nc=f_jhi_nc, beta=tempBeta, zeta=tempZeta)
end