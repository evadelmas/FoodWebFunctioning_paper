#=
Feb. 26th 2020
Eva Delmas

This script will generate the interaction matrices of the food webs needed for
the analysis. We want to have food webs representing communities with
- different levels of diversity (producer and consumer richness)
- different structures for the same level of diversity
We use the Allometric Diet Breadth Model (presented in Petchey et al., 2008) to
simulate realistic food webs. As we want to compare the effect of richness and
not ecosystem type, we will sample the body masses needed to initialize the
model  in a single community: the Benguela pelagic ecosystem. The data are in
the  data/benguela_bodysizes.txt file.
The functions implementing the ADBM model are in scripts/ADBM_model.jl
=#
using Distributed
addprocs(8)

@everywhere import Pkg
@everywhere Pkg.activate(pwd())
@everywhere Pkg.instantiate()
@everywhere include("scripts/ADBM_model.jl")
@everywhere import BioEnergeticFoodWebs.trophic_rank
using JLD2

P = collect(2:2:20) #plant diversity levels
C = collect(5:5:35) #consumers diversity levels
PC = [] #combine
for i in P ; for j in C ; append!(PC, [(i,j)]) ; end ; end
n_structures = 100
fw_todo = repeat(PC, n_structures) #we want to have n_structures food webs that have the same diversity but contain different species, organized possibly in different ways

matrices = pmap(x -> Main.adbm_model(x[1], x[2]), fw_todo) #generating the food webs
matrices_unique = unique(matrices) #removing duplicates
#NB: The names of the rows contain a number to identify the species and its bodysize,
#the name of the colums contain its metabolic status (false = vertebrate, true = invertebrate)
#as such, it can be retrieve if needed for initializing the simulations

#Check for disconnected species, food webs without producers and food webs with
#no more than 2 trophic levels (herbivores only).

tmp = Array.(matrices_unique)
for i in 1:length(tmp)
    A = tmp[i]
    is_connected = vec(sum(A, dims = 2) .!= 0) .| vec(sum(A, dims = 1) .!= 0)
    C = all(is_connected)
    is_producer = vec(sum(A, dims = 2) .== 0)
    P = any(is_producer)
    height = trophic_rank(A)
    H = any(height .> 2)
    !(C & P & H) ? deleteat!(matrices_unique, i) : Nothing
end

@save "data/adbm_generated_foodwebs.jld2" matrices_unique
