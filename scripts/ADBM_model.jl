#=
This scripts contains the functions needed to simulate food webs based on the
Allometric Diet Breadth Model (ADBM -- Petchey et al., 2008).
Author: Eva Delmas
Date: June, 2020
=#

using NamedArrays, DelimitedFiles
import StatsBase.sample
import BioEnergeticFoodWebs.trophic_rank

#=
Load and clean raw data from https://figshare.com/collections/BODY_SIZES_OF_CONSUMERS_AND_THEIR_RESOURCES/3298772
=#
using DelimitedFiles, DataFrames, CSVFiles
raw_data_all = readdlm("data/raw_data/bodysizes_2008.txt", '\t', header = true)
raw_data_data = raw_data_all[1]
header_raw_data = raw_data_all[2]
f(A,names) = DataFrame(Any[A[:,i] for i = 1:size(A,2)], map(Symbol,names))
raw_data = f(raw_data_data, vec(header_raw_data))
#restrict to Benguela
idx_benguela = raw_data[!,Symbol("Geographic location")] .== "Africa, Benguela ecosystem"
idx_nobact = raw_data[!,Symbol("Type of feeding interaction")] .!= "bacterivorous"
benguela_all = raw_data[idx_benguela .& idx_nobact,:]
resource_isvert = (benguela_all[!,Symbol("Metabolic category resource")] .== "ectotherm vertebrate") .| (benguela_all[!,Symbol("Metabolic category resource")] .== "endotherm vertebrate")
consumer_isvert = (benguela_all[!,Symbol("Metabolic category consumer")] .== "ectotherm vertebrate") .| (benguela_all[!,Symbol("Metabolic category consumer")] .== "endotherm vertebrate")
resource_isprod = (benguela_all[!,Symbol("Metabolic category resource")] .== "photo-autotroph")
benguela = DataFrame(
   interaction_type = benguela_all[!,Symbol("Type of feeding interaction")]
   , consumer_species = benguela_all[!,Symbol("Common name(s) consumer")]
   , consumer_metab = benguela_all[!,Symbol("Metabolic category consumer")]
   , consumer_mass = benguela_all[!,Symbol("Mean mass (g) consumer")]
   , resource_species = benguela_all[!,Symbol("Common name(s) resource")]
   , resource_metab = benguela_all[!,Symbol("Metabolic category resource")]
   , resource_mass = benguela_all[!,Symbol("Mean mass (g) resource")]
)
producers_bm = benguela[resource_isprod,:resource_mass]
consumers_bm = vcat(benguela[.!resource_isprod,:resource_mass], benguela[!,:consumer_mass])
isvert = vcat(resource_isvert[.!resource_isprod], consumer_isvert)
#CSVFiles.save("data/benguela_cleandata.csv", benguela)

function bodysize_empirical(P::Int64, C::Int64; producers_bm=producers_bm, consumers_bm=consumers_bm, is_vertebrate = isvert)
   #bsize_data = readdlm("data/benguela_bodysizes.txt") #open body size file
   bsize_out = zeros(P+C)
   isvert_out = falses(P+C)
   #autotrophs = rand(findall((bsize_data[!,:] .== "photo-autotroph") .& (bsize_data[:,5] .!= -999)), P)
   #consumers = rand(findall((bsize_data[:,2] .!= "photo-autotroph") .& (bsize_data[:,3] .!= -999)), C)
   consumers = rand(1:length(is_vertebrate),C)
   bsize_out[1:P] = rand(producers_bm, P)
   bsize_out[P+1:end] = consumers_bm[consumers]
   isvert_out[P+1:end] = is_vertebrate[consumers]
   return [bsize_out, isvert_out]
end

function bodysize_sample(P::Int64, C::Int64)
   bsize_out = zeros(P+C)
   distr_autotrophs = TruncatedNormal(0,0.1,1e-10,1)
   distr_heterotrophs = TruncatedNormal(0,50,1e-10,150)
   bsize_out[1:P] = rand(distr_autotrophs, P)
   bsize_out[P+1:end] = rand(distr_heterotrophs, C)
   return bsize_out
end

function adbm_internal(P::Int64, C::Int64; bm::Array{Float64,1} = Array{Float64,1}(undef, 0))
   if bm == Array{Float64,1}(undef, 0)
      emp_data = bodysize_empirical(P, C)
      bm = emp_data[1]
      isvert = emp_data[2]
   end
   @assert length(bm) == P+C
   S = P + C
   E = bm
   N = bm .^ -0.75
   A = 0.0189 * (bm .^ -0.465) * (bm .^ -0.491)'
   for i = 1:S #for each prey
      A[:,i] = A[:,i] .* N[i]
   end
   λ = A
   H = zeros(Float64,size(A))
   ratios = (bm ./ bm')'
   for i = 1:S , j = 1:S
      if ratios[j,i] < 0.401
         H[j,i] =  1 / (0.401 - ratios[j,i])
      else
         H[j,i] = Inf
      end
   end
   adbmMAT = zeros(Int64,size(A))
   for i in 1:S
      feeding = get_feeding_links(S, E, λ, H, bm, i)
      if length(feeding) > 0
         for j in feeding
            adbmMAT[i,j] = 1
         end
      end
   end
   bm_sorted = sortperm(bm)
   adbmMAT_s = NamedArray(adbmMAT[bm_sorted,bm_sorted])
   linenames = string.(1:1:S) .* "_" .* string.(bm[bm_sorted])
   colnames = string.(1:1:S) .* "_" .* string.(isvert[bm_sorted])
   setnames!(adbmMAT_s, linenames, 1)
   setnames!(adbmMAT_s, colnames, 2)
   return adbmMAT_s
end

function get_feeding_links(S::Int64,E::Vector{Float64}, λ::Matrix{Float64},
   H::Matrix{Float64},bm::Vector{Float64},j)

   profit = E ./ H[j,:]
   # Setting profit of species with zero biomass  to -1.0
   # This prevents them being included in the profitSort
   profit[bm .== 0.0] .= -1.0

   profs = sortperm(profit,rev = true)

   λSort = λ[j,profs]
   HSort = H[j,profs]
   ESort = E[profs]

   λH = cumsum(λSort .* HSort)
   Eλ = cumsum(ESort .* λSort)

   λH[isnan.(λH)] .= Inf
   Eλ[isnan.(Eλ)] .= Inf

   cumulativeProfit = Eλ ./ (1 .+ λH)

   if all(0 .== cumulativeProfit)
      feeding = []
   else
      feeding = profs[1:maximum(findall(cumulativeProfit .== maximum(cumulativeProfit)))]
   end

   #cumulativeProfit[end] = NaN
   #feeding = profs[(append!([true],cumulativeProfit[1:end-1] .< profitSort[2:end]))]
   return(feeding)
end

function make_one_adbm_attempt!(A, P::Int64, C::Int64; bodymasses::Array{Float64,1} = Array{Float64,1}(0))
   for i in eachindex(A)
      A[i] = 0
   end
   potential_links = adbm_internal(P, C, bm = bodymasses)
   link_pos_cart = findall(potential_links .== 1)
   link_pos = hcat((i->i[1]).(link_pos_cart), (i->i[2]).(link_pos_cart))
   link_pos = link_pos[link_pos[:,1] .> P,:]
   min_L = max(P, C)
   max_L = size(link_pos,1)
   setnames!(A, names(potential_links)[1], 1)
   setnames!(A, names(potential_links)[2], 2)
   if max_L >= min_L
      L = rand(min_L:max_L)
      new_links = link_pos[sample(1:size(link_pos,1), L, replace = false),:]
      for i in 1:size(new_links,1)
         x,y = new_links[i,:]
         A[x,y] = 1
      end
   end
end

function valid_adbm(A, P)
   return (sum(sum(A, dims = 2).==0) == P) & (minimum(sum(A, dims = 1).+sum(A, dims = 2)') > 0)
end

function adbm_model(P::Int64, C::Int64; bodymasses::Vector{Float64} = Array{Float64,1}(undef, 0))
   #println(string(P) * "_" * string(C))
   s = (P+C,P+C)
   A = NamedArray(zeros(Int64, s))
   make_one_adbm_attempt!(A, P, C, bodymasses=bodymasses)
   while !valid_adbm(A,P)
      make_one_adbm_attempt!(A, P, C, bodymasses=bodymasses)
   end
   return A
end

function get_bodysize(A, metab, Z_vert, Z_invert)
    p = sum(A,2) .== 0
    TR = trophic_rank(A)
    M = zeros(size(A,1))
    M[metab] = Z_vert .^ (TR[metab] .- 1)
    M[vcat((.!metab) .& (.!p)...)] = Z_invert .^ (TR[vcat((.!metab) .& (.!p)...)] .- 1)
    M[vcat(p...)] = 1.0
    return M
end
