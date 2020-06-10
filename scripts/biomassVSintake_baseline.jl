#=
March 7th 2020
Eva Delmas

This script runs the simulations necessary to estimate the baseline
for the biomass to intake relationship (needed for the figures).
This baseline is in fact the relationship between the total biomass
and the total intake of a prairie (plants only).
=#
using BioEnergeticFoodWeb
using DelimitedFiles

import StatsBase.sample
N = collect(1:2:10) #supply range
S = collect(1:4:20) #richness range
ks = collect(0.15:0.01:0.2)
comb = []
for n in N
    for s in S
        for k in ks
            append!(comb, [[n, s, k]])
        end
    end
end
function simulate_intake(c)
    S = Int(c[2])
    A = Int.(zeros(S,S))
    k1 = k2 = repeat([c[3]], S)
    N = c[1]
    p = model_parameters(A, productivity = :nutrients, K1 = k1, K2 = k2, supply = [N, N])
    out = simulate(p, rand(S), stop = 2000)
    t_total = length(out[:t])
    if t_total > 8000
        tkeep = Int(round(t_total / 3))
        intake = sum(nutrient_intake(out, last = tkeep, out_type = :mean))
        biomass = total_biomass(out, last = tkeep)
        return (intake, biomass)
    end
end

res = Array{Float64,2}(undef,length(comb),2)
for i in 1:length(comb)
    res[i,1], res[i,2] = simulate_intake(comb[i])
end

writedlm("simulations_outputs/baseline.csv", res)
