#=
Feb. 26th 2020
Eva Delmas

=#

#REQUIREMENTS ==================================================================

using Distributed
addprocs(8)
@everywhere import Pkg
@everywhere Pkg.activate(pwd())
@everywhere using JLD2
@everywhere include("scripts/utils.jl")
using DataFrames, CSV

#EXTRACTS RESULTS ==============================================================

out_files = readdir("simulations_outputs/")
out_files = "simulations_outputs/" .* out_files

@everywhere function extract_results(simfile)
    @load simfile sim
    id = (id = string(hash(sim[:p][:A])),)
    tk = Int(round(length(sim[:t]) / 3))

    B = biomass(sim, tkeep = tk)
    F = flux(sim, tkeep = tk)
    St = nkstructure(sim)
    Prop = properties(sim)
    Stab = stability(sim)
    M = motifs(sim)
    Motifs = NamedTuple{Tuple(keys(M))}(values(M))
    SR = species_roles(sim)
    S = size(sim[:p][:newA], 1)
    Div = (S = S
        , P = sum(SR.Producers)
        , H = sum(SR.Herbivores)
        , O = sum(SR.Omnivores)
        , C = sum(SR.Carnivores))

    out = merge(id, Div, B, F, St, Prop, Stab, Motifs)
end

results = pmap(extract_results, out_files)
results_df = DataFrame(results)
CSV.write("data/results.csv", results_df)
