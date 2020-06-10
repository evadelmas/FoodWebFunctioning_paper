#=
March 3rd 2020
Eva Delmas

This script contains all the functions needed for extracting measures of food-web
structure, properties, ets.
=#

using StatsBase, BioEnergeticFoodWebs, DelimitedFiles, PyCall

function species_roles(sim)
    #NB: here omnivores are species that eat both consumers (primary or secondary)
    #and producers.
    A = sim[:p][:newA]
    TL = trophic_rank(A)
    P = vec(sum(A, dims = 2) .== 0)
    eat_prod = vec(sum(A[:,vec(P)], dims = 2) .!= 0)
    eat_no_prod = vec(sum(A[:,vec(P)], dims = 2) .== 0)
    eat_no_nonprod = vec(sum(A[:,vec(.!P)], dims = 2) .== 0)
    eat_nonprod = vec(sum(A[:,vec(.!P)], dims = 2) .!= 0)
    H = vec(eat_prod .& eat_no_nonprod)
    O = vec(eat_prod .& eat_nonprod)
    C = vec(eat_no_prod .& eat_nonprod)
    #test
    @assert all(TL[H] .== 2.0)
    @assert all(TL[C] .> 2.0)
    @assert all(TL[O] .> 2.0)
    @assert all(TL[P] .== 1.0)
    return (Producers = P, Herbivores = H, Carnivores = C, Omnivores = O, TL = TL)
end

function trophic_structure(bcomp)
    dominant_bm = findmax(bcomp)[2]
    order_tl = sortperm(bcomp, rev = true)
    if (order_tl == [1,2,3]) .| (order_tl == [3,2,1])
        shape = "Pyramid"
    else
        shape = "Cascade"
    end
    return(dominant_bm, shape)
end

function biomass(sim; tkeep::Int64 = 6000)
    P, H, C, O, TL = species_roles(sim)
    A = sim[:p][:newA]
    nonextinct = trues(sim[:p][:S])
    nonextinct[sim[:p][:extinctions]] .= false
    b = population_biomass(sim, last = tkeep)[nonextinct]
    top_heaviness = (b ./ b') .* A
    th_herb = (b ./ b')[H .| O, P] .* A[H .| O, P]
    th_pred = (b ./ b')[C .| O, .!P] .* A[C .| O, .!P]
    TH = mean(top_heaviness[findall(A .!= 0)])
    TH_herb = mean(th_herb[findall(A[H .| O, P].!= 0)])
    TH_pred = mean(th_pred[findall(A[C .| O, .!P].!= 0)])
    B = sum(b)
    b_h = sum(b[H])
    b_p = sum(b[P])
    b_o = sum(b[O])
    b_c = sum(b[C])
    @assert isapprox(B, (b_h + b_p + b_o + b_c), atol = 1e-6)
    dom, shape = trophic_structure([b_p, b_h, b_o + b_c])
    b23 = sum(b[(TL .> 2) .& (TL .< 3)])
    b3 = sum(b[(TL .== 3)])
    b34 = sum(b[(TL .> 3) .& (TL .< 4)])
    b4 = sum(b[(TL .== 4)])
    b45 = sum(b[(TL .> 4) .& (TL .< 5)])
    b5 = sum(b[(TL .== 5)])
    b56 = sum(b[(TL .> 5)])
    @assert isapprox(B, (b_p + b_h + b23 + b3 + b34 + b4 + b45 + b5 + b56), atol = 1e-6)
    return(b = B, bP = b_p, bH = b_h, bO = b_o, bC = b_c,
        shape = shape, dom = dom, TH = TH, THherb = TH_herb, THpred = TH_pred,
        b23 = b23, b3 = b3, b34 = b34, b4 = b4, b45 = b45, b5 = b5, b56 = b56)
end

function flux_matrix(sim; tkeep::Int64 = 6000)
    bm = population_biomass(sim, last = tkeep)
    parameters = sim[:p]

    # Total available biomass
    bm_matrix = zeros(eltype(bm), (length(bm), length(bm)))
    rewire = (parameters[:rewire_method] == :ADBM) | (parameters[:rewire_method] == :Gilljam) | (parameters[:rewire_method] == :DS)
    costMat = rewire ? parameters[:costMat] : nothing
    BioEnergeticFoodWebs.fill_bm_matrix!(bm_matrix, bm, parameters[:w], parameters[:A], parameters[:h]; rewire=rewire, costMat=costMat)

    # Available food
    F = zeros(eltype(bm), (length(bm), length(bm)))
    BioEnergeticFoodWebs.fill_F_matrix!(F, bm_matrix, bm, parameters[:Γh], parameters[:c])

    # XYB matrix
    xyb = zeros(eltype(bm), length(bm))
    BioEnergeticFoodWebs.fill_xyb_matrix!(xyb, bm, parameters[:x], parameters[:y])

    BioEnergeticFoodWebs.update_F_matrix!(F, xyb)
    in_mat = deepcopy(F)
    BioEnergeticFoodWebs.get_trophic_loss!(F, parameters[:efficiency])
    out_mat = deepcopy(F)
    return in_mat, out_mat
end

function flux(sim; tkeep::Int64 = 6000)
    P, H, C, O, TL = species_roles(sim)
    nonextinct = trues(sim[:p][:S])
    nonextinct[sim[:p][:extinctions]] .= false
    intake_all = vec(nutrient_intake(sim; last = tkeep, out_type = :mean))
    flux_in, flux_out = flux_matrix(sim, tkeep = tkeep)
    flux_in = flux_in[nonextinct, nonextinct]
    flux_out = flux_out[nonextinct, nonextinct]
    flux_in_ph = []
    flux_in_po = []
    flux_in_ho = []
    flux_in_hc = []
    flux_in_oo = []
    flux_in_cc = []
    flux_in_oc = []
    flux_in_co = []
    flux_out_ph  = []
    flux_out_po  = []
    flux_out_ho  = []
    flux_out_hc  = []
    flux_out_oo  = []
    flux_out_cc  = []
    flux_out_oc  = []
    flux_out_co  = []
    for c in 1:length(TL)
        for r in 1:length(TL)
            if (r ∈ findall(P)) & (c ∈ findall(H))
                push!(flux_in_ph, flux_in[c,r])
                push!(flux_out_ph, flux_out[c,r])
            elseif (r ∈ findall(P)) & (c ∈ findall(O))
                push!(flux_in_po, flux_in[c,r])
                push!(flux_out_po, flux_out[c,r])
            elseif (r ∈ findall(H)) & (c ∈ findall(O))
                push!(flux_in_ho, flux_in[c,r])
                push!(flux_out_ho, flux_out[c,r])
            elseif (r ∈ findall(H)) & (c ∈ findall(C))
                push!(flux_in_hc, flux_in[c,r])
                push!(flux_out_hc, flux_out[c,r])
            elseif (r ∈ findall(O)) & (c ∈ findall(O))
                push!(flux_in_oo, flux_in[c,r])
                push!(flux_out_oo, flux_out[c,r])
            elseif (r ∈ findall(C)) & (c ∈ findall(O))
                push!(flux_in_co, flux_in[c,r])
                push!(flux_out_co, flux_out[c,r])
            elseif (r ∈ findall(C)) & (c ∈ findall(C))
                push!(flux_in_cc, flux_in[c,r])
                push!(flux_out_cc, flux_out[c,r])
            elseif (r ∈ findall(O)) & (c ∈ findall(C))
                push!(flux_in_oc, flux_in[c,r])
                push!(flux_out_oc, flux_out[c,r])
            end
        end
    end
    m_all = metabolism(sim, out_type = :mean)[nonextinct]
    cons = sum(flux_in)
    return (cons = cons
        , Pin = length(intake_all) != 0 ? sum(intake_all) : 0.0
        , PHin = length(flux_in_ph) != 0 ? sum(flux_in_ph) : 0.0
        , PHout = length(flux_out_ph) != 0 ? sum(flux_out_ph) : 0.0
        , POin = length(flux_in_po) != 0 ? sum(flux_in_po) : 0.0
        , POout = length(flux_out_po) != 0 ? sum(flux_out_po) : 0.0
        , HOin = length(flux_in_ho) != 0 ? sum(flux_in_ho) : 0.0
        , HOout = length(flux_out_ho) != 0 ? sum(flux_out_ho) : 0.0
        , HCin = length(flux_in_hc) != 0 ? sum(flux_in_hc) : 0.0
        , HCout = length(flux_out_hc) != 0 ? sum(flux_out_hc) : 0.0
        , OOin = length(flux_in_oo) != 0 ? sum(flux_in_oo) : 0.0
        , OOout = length(flux_out_oo) != 0 ? sum(flux_out_oo) : 0.0
        , CCin = length(flux_in_cc) != 0 ? sum(flux_in_cc) : 0.0
        , CCout = length(flux_out_cc) != 0 ? sum(flux_out_cc) : 0.0
        , OCin = length(flux_in_oc) != 0 ? sum(flux_in_oc) : 0.0
        , OCout = length(flux_out_oc) != 0 ? sum(flux_out_oc) : 0.0
        , COin = length(flux_in_co) != 0 ? sum(flux_in_co) : 0.0
        , COout = length(flux_out_co) != 0 ? sum(flux_out_co) : 0.0
        , mP = length(m_all[P]) != 0 ? sum(m_all[P]) : 0.0
        , mH = length(m_all[H]) != 0 ? sum(m_all[H]) : 0.0
        , mO = length(m_all[O]) != 0 ? sum(m_all[O]) : 0.0
        , mC = length(m_all[C]) != 0 ? sum(m_all[C]) : 0.0)
end

function nkstructure(sim)
    A = sim[:p][:newA]
    P, H, C, O, TL = species_roles(sim)

    #connectance
    L = sum(A)

    #link between levels
    PH = sum(A[H,P])
    PO = sum(A[O,P])
    HO = sum(A[O,H])
    OO = sum(A[O,O])
    OC = sum(A[C,O])
    CO = sum(A[O,C])
    CC = sum(A[C,C])
    HC = sum(A[C,H])

    #height
    height = maximum(TL)
    height_mean = mean(TL)

    return (L = L
        , lPH = PH
        , lPO = PO
        , lHO = HO
        , lOO = OO
        , lOC = OC
        , lCO = CO
        , lCC = CC
        , lHC = HC
        , Height = height
        , Height_mean = height_mean)
end

function properties(sim)
    P, H, C, O, TL = species_roles(sim)
    p = sim[:p]
    nonextinct = trues(p[:S])
    nonextinct[p[:extinctions]] .= false
    A = p[:newA]

    #bodymass
    bodymass = p[:bodymass][nonextinct]
    ratio_bm = (bodymass ./ bodymass') .* A
    bm_herb = (bodymass ./ bodymass')[H .| O, P] .* A[H .| O, P]
    bm_pred = (bodymass ./ bodymass')[C .| O, .!P] .* A[C .| O, .!P]
    Z = mean(ratio_bm[findall(A .!= 0)])
    Z_herb = mean(bm_herb[findall(A[H .| O, P].!= 0)])
    Z_pred = mean(bm_pred[findall(A[C .| O, .!P].!= 0)])
    bPmean = sum(P) != 0 ? mean(bodymass[P]) : NaN
    bPmax = sum(P) != 0 ? maximum(bodymass[P]) : NaN
    bHmean = sum(H) != 0 ? mean(bodymass[H]) : NaN
    bHmax = sum(H) != 0 ? maximum(bodymass[H]) : NaN
    bOmean = sum(O) != 0 ? mean(bodymass[O]) : NaN
    bOmax = sum(O) != 0 ? maximum(bodymass[O]) : NaN
    bCmean = sum(C) != 0 ? mean(bodymass[C]) : NaN
    bCmax = sum(C) != 0 ? maximum(bodymass[C]) : NaN

    #metabolic class
    p_vertebrate = sum(p[:vertebrates][nonextinct]) / (sum(nonextinct) - sum(P))

    #producers efficiency
    k1 = p[:K1][nonextinct][P]
    k2 = p[:K2][nonextinct][P]
    k1min = minimum(k1)
    k2min = minimum(k2)
    k1mean = mean(k1)
    k2mean = mean(k2)
    k1max = maximum(k1)
    k2max = maximum(k2)

    return (Z = Z
        , Zh = Z_herb
        , Zp = Z_pred
        , bPmean = bPmean
        , bPmax = bPmax
        , bHmean = bHmean
        , bHmax = bHmax
        , bOmean = bOmean
        , bOmax = bOmax
        , bCmean = bCmean
        , bCmax = bCmax
        , p_vertebrate = p_vertebrate
        , k1min = k1min
        , k2min = k2min
        , k1mean = k1mean
        , k2mean = k2mean
        , k1max = k1max
        , k2max = k2max)
end

function motifs(sim; tkeep::Int64 = 6000)
    A = sim[:p][:newA]
    function matrix_to_list(A)
        idx = findall(A .== 1)
        from_sp = (i->i[1]).(idx)
        to_sp = (i->i[2]).(idx)
        return hcat(from_sp, to_sp)
    end
    webname = "tmpdata/" * string(hash(A)) * ".txt"
    writedlm(webname, matrix_to_list(A))
    py"""
    import pymfinder as pf

    def mc(x):
        return str(pf.motif_structure(x, stoufferIDs = True, allmotifs = True))
    """
    local motifs_counts
    try
        output = py"mc"(webname)
        out_split = split(output, "\n")
        out_split = out_split[length.(out_split) .!= 0]
        out_split_all = split.(out_split, " ")
        motifs_counts = Dict(:LFC => 0, :EC => 0, :AC => 0, :OMN => 0, :LOOP => 0, :TOT => 0)
        motifs_counts[:TOT] = sum(map(x -> parse(Int64, x[2]), out_split_all[2:end]))
        for m in out_split_all[2:end]
            if m[1] == "S1"
                motifs_counts[:LFC] = parse(Int64, m[2])
            elseif m[1] == "S2"
                motifs_counts[:OMN] = parse(Int64, m[2])
            elseif m[1] == "S3"
                motifs_counts[:LOOP] = parse(Int64, m[2])
            elseif m[1] == "S4"
                motifs_counts[:EC] = parse(Int64, m[2])
            elseif m[1] == "S5"
                motifs_counts[:AC] = parse(Int64, m[2])
            end
        end
    catch
        motifs_counts = Dict(:LFC => NaN, :EC => NaN, :AC => NaN, :OMN => NaN, :LOOP => NaN, :TOT => NaN)
    end

    return motifs_counts
end

function stability(sim; tkeep::Int64 = 6000)
    cv = population_stability(sim, last = tkeep)
    evenness = foodweb_evenness(sim, last = tkeep)
    persistence = size(sim[:p][:newA], 1) / size(sim[:p][:A], 1)

    return (cv = cv, evenness = evenness, persistence = persistence)
end
