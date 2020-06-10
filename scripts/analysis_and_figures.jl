fw#=
Feb. 26th 2020
Eva Delmas

This script generates the figures for the manuscript.
All figures are saved to manuscript/figures.
=#

#REQUIREMENTS ==================================================================

import Pkg
Pkg.activate(pwd())
Pkg.instantiate()
using Plots, Statistics, StatsPlots, DataFrames, CSVFiles, FileIO, DataFrames, DelimitedFiles

fw_data = load("data/results.csv") |> DataFrame

#get one sim for each web
id_unique = unique(fw_data[!,:id])
fw_idx = []
for i in id_unique[2:end]
    tmp = findall(fw_data[!,:id] .== i)
    append!(fw_idx, rand(tmp))
end
fw = fw_data[fw_idx,:]

#fonts
fntlg = Plots.font(16)
fntlg_r = deepcopy(fntlg)
fntlg_r.rotation = 90

#colors
cols = Dict(:green => RGB(0/255,158/255,115/255),
            :black => RGB(0,0,0),
            :blue => RGB(0/255,114/255,178/255),
            :vermillon => RGB(213/255,94/255,0/255),
            :yellow => RGB(240/256, 228/256, 66/256),
            :skyblue => RGB(86/256, 180/256, 233/256),
            :purple => RGB(204/256, 121/256, 167/256),
            :orange => RGB(230/256, 159/256, 0/256))
cols_light = Dict(:green => RGB(204/255,255/255,241/255),
                  :black => RGB(204/255,204/255,204/255),
                  :blue => RGB(204/255,236/255,255/255),
                  :orange => RGB(255/255,226/255,204/255))

#=
    FIGURES ====================================================================
=#

dom = fw[!,:dom]

#Biomass to intake
import DelimitedFiles.readdlm
import DelimitedFiles.writedlm
res = readdlm("data/baseline.csv")
using GLM
res_df = DataFrame(x = res[:,1], y = res[:,2])
ols = lm(@formula(y ~ x), res_df)
coefficients = coef(ols)
keep_coef = (coef(ols) .- stderror(ols) .<= 0) .& (coef(ols) .+ stderror(ols) .>= 0)
coefficients[keep_coef] .= 0.0

baseline(x) = coefficients[1] + coefficients[2] * x

intake = fw[!,:Pin]
base_biomass = baseline.(intake)
is_below = round(sum(fw[!,:b] .< base_biomass) / length(intake) *100, digits = 1)
is_above = round(sum(fw[!,:b] .> base_biomass) / length(intake) *100, digits = 1)

#write the number of food web in each case in txt files
writedlm("manuscript/figures/nfw.txt", size(fw,1))
writedlm("manuscript/figures/below.txt", is_below)
writedlm("manuscript/figures/above.txt", is_above)
writedlm("manuscript/figures/nbh.txt", round(sum(dom .== 1) / size(fw,1) *100, digits = 1))
writedlm("manuscript/figures/nth.txt", round(sum(dom .== 3) / size(fw,1) *100, digits = 1))
writedlm("manuscript/figures/nmh.txt", round(sum(dom .== 2) / size(fw,1) *100, digits = 1))
writedlm("manuscript/figures/nbhc.txt", round(sum((dom .== 1) .& (fw[!,:shape] .== "Cascade")) / sum(dom .== 1) *100, digits = 1))
writedlm("manuscript/figures/nbhp.txt", round(sum((dom .== 1) .& (fw[!,:shape] .== "Pyramid")) / sum(dom .== 1) *100, digits = 1))
writedlm("manuscript/figures/nthc.txt", round(sum((dom .== 3) .& (fw[!,:shape] .== "Cascade")) / sum(dom .== 3) *100, digits = 1))

# FIGURE 4 - Biomass to intake
p = []
ltr = ["A", "B", "C", "D", "E"]

global i = 0
for d in [3,2,1]
    for shape in ["Cascade", "Pyramid"]
        toplot = (fw[!,:shape] .== [shape]) .& (dom .== d)
        N = sum(toplot)
        if sum(toplot) != 0
            global i = i .+ 1
            lt = ltr[i]
            tmp = fw[vec(toplot),:]
            #tmpdf = DataFrame(x = round.(tmp[!,:Pin], digits = 2), y = tmp[!,:b])
            #cons_mean = aggregate(tmpdf, :x, [mean, std])
            lmar = shape == "Cascade" ? 15Plots.mm : 7Plots.mm
            pi = scatter(tmp[!,:Pin], tmp[!,:b]
                #, mc = :grey
                , msw = 0.5, ma = 0.3, msa = 0.1
                , c = :BrBG
                , mz = log10.(tmp[!,:Z]), clim = (-4.7,4.7)
                , ylims = (0,14), xlims = (0.0,1.01)
                , leg = false
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                #, yerr = assim_mean[!,:y_std]
                )
            xticks!([0.0:0.2:1;])
            yticks!([0:2:12;])
            annotate!(0.1, 13, text("$lt", fntlg))
            tmp_i = collect(0:0.2:1.0)
            # savefig("figure/regime_shape$i$j" * ".png")
            plot!(tmp_i, baseline, linecolor = :black, linestyle = :dash, lw = 0.5, labels = "")
            d == 1 ? xaxis!("total intake") : Nothing
            shape == "Cascade" ? yaxis!("total biomass") : Nothing
            if d == 3
                if shape == "Cascade"
                    title!("CASCADE", titlefont = fntlg)
                else
                    title!("PYRAMID", titlefont = fntlg)
                end
            end
            push!(p, pi)
            b = DataFrame(P = tmp[!,:bP], H = tmp[!,:bH], C = tmp[!,:bO] .+ tmp[!,:bC])
            pj = plot(legend=false,grid=false,foreground_color_border=:white)
            for c in 1:size(b,2)
                bm = mean(b[:,c])
                plgn = Shape([-bm/2,bm/2,bm/2,-bm/2], [0+c,0+c,1+c,1+c])
                plot!(plgn, fill = (0, 0.3, :grey), linecolor = :grey, left_margin = lmar, xlims = (-4,4), ylims = (1,5))
            end
            annotate!(0, 4.5, text("N = $N", fntlg))
            xaxis!(false)
            yticks!([1.5:1:3.5;], ["P", "C1", "C2"])
            if shape == "Cascade"
                if d == 1
                    annotate!(-8, 2.5, text("BOTTOM HEAVY", fntlg_r))
                elseif d == 2
                    annotate!(-8, 2.5, text("MIDDLE HEAVY", fntlg_r))
                else
                    annotate!(-8, 2.5, text("TOP HEAVY", fntlg_r))
                end
            end
            push!(p, pj)
        else
            pi = scatter([NaN],[NaN],grid=false
                , xaxis = false, yaxis = false
                , xlims = (0,1), ylims = (0,1)
                , mc = :grey, ma = 0.4, msw = 0.5
                , legend = :topright, label = "Food-web biomass"
                , legendfont = Plots.font(12)
                , foreground_color_legend = nothing)
            plot!([NaN], [NaN]
                , lc = :black, linestyle = :dash
                , lw = 0.5, labels = "No consumers")
            pj = plot([NaN],[NaN],grid=false
                , mz = fw[!,:Z], clim = (-4.7,4.7)
                , xaxis = false, yaxis = false
                , xlims = (0,1), ylims = (0,1)
                , c = :BrBG
                , legend = :best, label = ""
                , legendfont = Plots.font(14)
                , colorbar_title = "log10.(Z)"
                , foreground_color_legend = nothing)
            push!(p, pi)
            push!(p, pj)
        end
    end
end
l = grid(3,4, widths = [.15,.35,.35,.15])
porder = p[vcat(2,1,3,4,6,5,7,8,10,9,11,12)]
plot(porder..., layout = l, size = (1000,970))
savefig("manuscript/figures/biomass-to-intake.png")

#FIGURE 5 - BEF, flux
p = []
ltr1 = ["A", "B", "C", "D", "E"]
global i = 0
for d in [3,2,1]
    for shape in ["Cascade", "Pyramid"]
        toplot = (fw[!,:shape] .== [shape]) .& (fw[!,:dom] .== d)
        N = sum(toplot)
        if sum(toplot) != 0
            #get letters for label
            global i = i .+ 1
            lti = ltr1[i]
            #data subset
            tmp = fw[vec(toplot),:]
            #calculate flux
            herb = tmp[!,:PHin] .+ tmp[!,:POin]
            cons = tmp[!,:HOin] .+ tmp[!,:HCin]
            igp = tmp[!,:OOin] .+ tmp[!,:CCin] .+ tmp[!,:OCin] .+ tmp[!,:COin]
            tot = herb .+ cons .+ igp
            #calculate the mean and std
            tmpdf = DataFrame(A = tmp[!,:S] .- tmp[!,:P], P = tmp[!,:P], b = tmp[!,:b], h = herb, c = cons, o = igp, t = tot, i = tmp[!,:Pin])
            cons_mean_A = aggregate(tmpdf, :A, [mean, std], skipmissing = true)
            cons_mean_P = aggregate(tmpdf, :P, [mean, std], skipmissing = true)
            deleterows!(cons_mean_A, findall(isnan.(cons_mean_A[!,:h_std])))
            deleterows!(cons_mean_P, findall(isnan.(cons_mean_P[!,:h_std])))
            sort!(cons_mean_A, :A)
            sort!(cons_mean_P, :P)

            #plot flux
            pi = plot(cons_mean_A[!,:A], cons_mean_A[!,:t_mean]
                , lc = cols[:purple]
                , ylims = (0,1.1), xlims = (0,35)
                , leg = false
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                , ribbon = cons_mean_A[!,:t_std]
                , fc = cols[:purple], fa = 0.2
                )
            annotate!(5, 1.05, text("$lti", 14))
            shape == "Pyramid" ? yticks!([0:0.25:1;], ["", "", ""]) : yticks!([0:0.25:1;], string.([0:0.25:1;]))
            shape == "Cascade" ? yaxis!("Intake") : Nothing
            d == 1 ? xaxis!("Animal richness") : Nothing
            d != 1 ? xticks!([0:5:35;], ["", "", ""]) : xticks!([0:5:35;], string.([0:5:35;]))
            plot!(cons_mean_A[!,:A], cons_mean_A[!,:i_mean], lc = cols[:green], ribbon = cons_mean_A[!,:i_std], fc = cols[:green], fa = 0.2)
            plot!(cons_mean_A[!,:A], cons_mean_A[!,:h_mean], lc = cols[:blue], ribbon = cons_mean_A[!,:h_std], fc = cols[:blue], fa = 0.2)
            plot!(cons_mean_A[!,:A], cons_mean_A[!,:c_mean], lc = cols[:vermillon], ribbon = cons_mean_A[!,:c_std], fc = cols[:vermillon], fa = 0.2)
            plot!(cons_mean_A[!,:A], cons_mean_A[!,:o_mean], lc = cols[:orange], ribbon = cons_mean_A[!,:o_std], fc = cols[:orange], fa = 0.2)

            #title
            if d == 3
                if shape == "Cascade"
                    title!("CASCADE", titlefont = fntlg)
                else
                    title!("PYRAMID", titlefont = fntlg)
                end
            end
            push!(p, pi)

            #shape
            b = DataFrame(P = tmp[!,:bP], H = tmp[!,:bH], C = tmp[!,:bO] .+ tmp[!,:bC])
            pk = plot(legend=false,grid=false,foreground_color_border=:white)
            lmar = shape == "Cascade" ? 15Plots.mm : 7Plots.mm
            for c in 1:size(b,2)
                bm = mean(b[:,c])
                plgn = Shape([-bm/2,bm/2,bm/2,-bm/2], [0+c,0+c,1+c,1+c])
                plot!(plgn, fill = (0, 0.3, :grey), linecolor = :grey, left_margin = lmar, xlims = (-4,4), ylims = (1,5))
            end
            annotate!(0, 4.5, text("N = $N", fntlg))
            xaxis!(false)
            yticks!([1.5:1:3.5;], ["P", "C1", "C2"], ytickfontsize = 12)

            #title2
            if shape == "Cascade"
                if d == 1
                    annotate!(-9, 2.5, text("BOTTOM  HEAVY", fntlg_r))
                elseif d == 2
                    annotate!(-9, 2.5, text("MIDDLE  HEAVY", fntlg_r))
                else
                    annotate!(-9, 2.5, text("TOP  HEAVY", fntlg_r))
                end
            end
            push!(p, pk)
        else

            #legend
            pk = plot([NaN NaN NaN NaN NaN],[NaN NaN NaN NaN NaN],grid=false
                , xaxis = false, yaxis = false
                , xlims = (0,1), ylims = (0,1)
                #, ls = [:solid :solid :dash :dash :dash :solid]
                , lc = [cols[:green] cols[:purple]  cols[:blue] cols[:vermillon] cols[:orange]], lw = 2
                , lab = ["Total plant intake" "Total animal intake" "Herbivory" "Carnivory (non IGP)" "Intraguild predation (IGP)"]
                , legend = :best
                , legendfont = fntlg
                , foreground_color_legend = nothing)
            pi = plot(legend=false,grid=false,foreground_color_subplot=:white, left_margin = 1Plots.mm)
            push!(p, pi)
            push!(p, pk)
        end
    end
end

l = @layout [
    a{0.15w} grid(1,2) c{0.15w}
    d{0.15w} grid(1,2) f{0.15w}
    g{0.15w} grid(1,2) h{0.15w}
]
p_order = [2,1,3,4,6,5,7,8,10,9,11,12]
plot(p[p_order]..., layout = l, size = (900,900))
savefig("manuscript/figures/bef_flux.png")

#FIGURE 6 - Contours BEF biomass
p = []
ltr1 = ["A", "B", "C", "D", "E"]
global i = 0
for d in [3,2,1]
    for shape in ["Cascade", "Pyramid"]
        toplot = (fw[!,:shape] .== [shape]) .& (dom .== d)
        #toplot = (dom .== d)
        N = sum(toplot)
        if sum(toplot) != 0
            #get letters for label
            global i = i .+ 1
            lti = ltr1[i]
            #data subset
            tmp = fw[vec(toplot),:]

            tmpdf = DataFrame(P = round.(tmp[!,:P] ./ 10) .* 10
                , A = round.((tmp[!,:S] .- tmp[!,:P]) ./ 5) .* 5
                , B = tmp[!,:b])
            b_mean = aggregate(tmpdf, [:A, :P], mean, skipmissing = true)

            x = b_mean[!,:P]
            y = b_mean[!,:A]
            z = b_mean[!,:B_mean]
            x_lvl = sort(unique(x))
            y_lvl = sort(unique(y))
            z_2d = zeros(length(y_lvl), length(x_lvl))
            for m in 1:length(y_lvl)
                for l in 1:length(x_lvl)
                    xval = x_lvl[l]:x_lvl[l] + 5
                    yval = y_lvl[m]:y_lvl[m] + 10
                    tmpz = [k ∈ xval for k in x] .& [k ∈ yval for k in y]
                    z_val = z[tmpz]
                    z_2d[m,l] = mean(z_val)
                end
            end
            thisp = contour(x_lvl, y_lvl, z_2d
                , clim = (0,10), lw = 1.5
                , ylims = (0,41), xlims = (0,21)
                , leg = false
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                )
            scatter!(tmp[!,:P], tmp[!,:S] .- tmp[!,:P], mz = tmp[!,:b]
                , clim = (0,10), msw = 0, ma = 0.3, ms = 3
                )
            annotate!(2, 38, text(lti, fntlg))
            d == 1 ? xaxis!("Producer richness") : Nothing
            shape == "Cascade" ? yaxis!("Consumer richness") : Nothing
            #annotate!(2,38, text(lti, 13))
            xticks!([0:5:20;], string.([0:5:20;]))
            yticks!([0:10:40;], string.([0:10:40;]))
            push!(p, thisp)

            if d == 3
                if shape == "Cascade"
                    title!("CASCADE", titlefont = fntlg)
                else
                    title!("PYRAMID", titlefont = fntlg)
                end
            end

            b = DataFrame(P = tmp[!,:bP], H = tmp[!,:bH], C = tmp[!,:bO] .+ tmp[!,:bC])
            pj = plot(legend=false,grid=false,foreground_color_border=:white)
            lmar = shape == "Cascade" ? 15Plots.mm : 7Plots.mm
            for c in 1:size(b,2)
                bm = mean(b[:,c])
                plgn = Shape([-bm/2,bm/2,bm/2,-bm/2], [0+c,0+c,1+c,1+c])
                plot!(plgn, fill = (0, 0.3, :grey), linecolor = :grey, left_margin = lmar, xlims = (-4,4), ylims = (1,5))
            end
            annotate!(0, 4.5, text("N = $N", fntlg))
            xaxis!(false)
            yticks!([1.5:1:3.5;], ["P", "C1", "C2"])
            if shape == "Cascade"
                if d == 1
                    annotate!(-8, 2.5, text("BOTTOM HEAVY", fntlg_r))
                elseif d == 2
                    annotate!(-8, 2.5, text("MIDDLE HEAVY", fntlg_r))
                else
                    annotate!(-8, 2.5, text("TOP HEAVY", fntlg_r))
                end
            end
            push!(p, pj)

        else

            #legend
            pj = plot([NaN],[NaN],grid=false
                , mz = fw[!,:Z], clim = (0,14)
                , xaxis = false, yaxis = false
                , xlims = (0,1), ylims = (0,1)
                , legend = :best, label = ""
                #, legendfont = fntlg
                , legendfontsize = 12
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                , colorbar_title = "Total biomass"
                , foreground_color_legend = nothing)
            push!(p, pj)
            pempty = plot([NaN],[NaN],grid=false,leg=false,axis=false)
            push!(p, pempty)
        end
    end
end
l = grid(3,4, widths = [.15,.35,.35,.15])
porder = p[vcat(2,1,3,4,6,5,7,8,10,9,11,12)]
plot(porder..., layout = l, size = (1000,970), bottom_margin = 15Plots.mm)
savefig("manuscript/figures/contour_bef.png")

#FIGURE S - contour biomass ratio
p = []
ltr1 = ["A", "B", "C", "D", "E"]
global i = 0
for d in [3,2,1]
    for shape in ["Cascade", "Pyramid"]
        toplot = (fw[!,:shape] .== [shape]) .& (dom .== d)
        #toplot = (dom .== d)
        N = sum(toplot)
        if sum(toplot) != 0
            #get letters for label
            global i = i .+ 1
            lti = ltr1[i]
            #data subset
            tmp = fw[vec(toplot),:]

            tmpdf = DataFrame(P = round.(tmp[!,:P] ./ 10) .* 10
                , A = round.((tmp[!,:S] .- tmp[!,:P]) ./ 5) .* 5
                , B = tmp[!,:bP] ./ tmp[!,:b])
            b_mean = aggregate(tmpdf, [:A, :P], mean, skipmissing = true)

            x = b_mean[!,:P]
            y = b_mean[!,:A]
            z = b_mean[!,:B_mean]
            x_lvl = sort(unique(x))
            y_lvl = sort(unique(y))
            z_2d = zeros(length(y_lvl), length(x_lvl))
            for m in 1:length(y_lvl)
                for l in 1:length(x_lvl)
                    xval = x_lvl[l]:x_lvl[l] + 5
                    yval = y_lvl[m]:y_lvl[m] + 10
                    tmpz = [k ∈ xval for k in x] .& [k ∈ yval for k in y]
                    z_val = z[tmpz]
                    z_2d[m,l] = mean(z_val)
                end
            end
            thisp = contour(x_lvl, y_lvl, z_2d
                , c = cgrad(:algae)
                , clim = (0.2,1), lw = 1.5
                , ylims = (0,41), xlims = (0,21)
                , leg = false
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                )
            scatter!(tmp[!,:P], tmp[!,:S] .- tmp[!,:P], mz = tmp[!,:bP] ./ tmp[!,:b]
                ,c = cgrad(:algae), clim = (0.2,1), msw = 0, ma = 0.3, ms = 3
                )
            annotate!(2, 38, text(lti, fntlg))
            d == 1 ? xaxis!("Producer richness") : Nothing
            shape == "Cascade" ? yaxis!("Consumer richness") : Nothing
            #annotate!(2,38, text(lti, 13))
            xticks!([0:5:20;], string.([0:5:20;]))
            yticks!([0:10:40;], string.([0:10:40;]))
            push!(p, thisp)

            if d == 3
                if shape == "Cascade"
                    title!("CASCADE", titlefont = fntlg)
                else
                    title!("PYRAMID", titlefont = fntlg)
                end
            end

            b = DataFrame(P = tmp[!,:bP], H = tmp[!,:bH], C = tmp[!,:bO] .+ tmp[!,:bC])
            pj = plot(legend=false,grid=false,foreground_color_border=:white)
            lmar = shape == "Cascade" ? 15Plots.mm : 7Plots.mm
            for c in 1:size(b,2)
                bm = mean(b[:,c])
                plgn = Shape([-bm/2,bm/2,bm/2,-bm/2], [0+c,0+c,1+c,1+c])
                plot!(plgn, fill = (0, 0.3, :grey), linecolor = :grey, left_margin = lmar, xlims = (-4,4), ylims = (1,5))
            end
            annotate!(0, 4.5, text("N = $N", fntlg))
            xaxis!(false)
            yticks!([1.5:1:3.5;], ["P", "C1", "C2"])
            if shape == "Cascade"
                if d == 1
                    annotate!(-8, 2.5, text("BOTTOM HEAVY", fntlg_r))
                elseif d == 2
                    annotate!(-8, 2.5, text("MIDDLE HEAVY", fntlg_r))
                else
                    annotate!(-8, 2.5, text("TOP HEAVY", fntlg_r))
                end
            end
            push!(p, pj)

        else

            #legend
            pj = plot([NaN],[NaN],grid=false
                , mz = fw[!,:Z], clim = (0.2,1)
                , c = cgrad(:algae)
                , xaxis = false, yaxis = false
                , xlims = (0,1), ylims = (0,1)
                , legend = :best, label = ""
                #, legendfont = fntlg
                , legendfontsize = 12
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                , colorbar_title = "b(prod) / b"
                , foreground_color_legend = nothing)
            push!(p, pj)
            pempty = plot([NaN],[NaN],grid=false,leg=false,axis=false)
            push!(p, pempty)
        end
    end
end
l = grid(3,4, widths = [.15,.35,.35,.15])
porder = p[vcat(2,1,3,4,6,5,7,8,10,9,11,12)]
plot(porder..., layout = l, size = (1000,970), bottom_margin = 15Plots.mm)
savefig("appendices/figures/contour_ratio.png")

#FIGURE S - contour productivity
p = []
ltr1 = ["A", "B", "C", "D", "E"]
global i = 0
for d in [3,2,1]
    for shape in ["Cascade", "Pyramid"]
        toplot = (fw[!,:shape] .== [shape]) .& (dom .== d)
        #toplot = (dom .== d)
        N = sum(toplot)
        if sum(toplot) != 0
            #get letters for label
            global i = i .+ 1
            lti = ltr1[i]
            #data subset
            tmp = fw[vec(toplot),:]

            tmpdf = DataFrame(P = round.(tmp[!,:P] ./ 10) .* 10
                , A = round.((tmp[!,:S] .- tmp[!,:P]) ./ 5) .* 5
                , B = tmp[!,:Pin])
            b_mean = aggregate(tmpdf, [:A, :P], mean, skipmissing = true)

            x = b_mean[!,:P]
            y = b_mean[!,:A]
            z = b_mean[!,:B_mean]
            x_lvl = sort(unique(x))
            y_lvl = sort(unique(y))
            z_2d = zeros(length(y_lvl), length(x_lvl))
            for m in 1:length(y_lvl)
                for l in 1:length(x_lvl)
                    xval = x_lvl[l]:x_lvl[l] + 5
                    yval = y_lvl[m]:y_lvl[m] + 10
                    tmpz = [k ∈ xval for k in x] .& [k ∈ yval for k in y]
                    z_val = z[tmpz]
                    z_2d[m,l] = mean(z_val)
                end
            end
            thisp = contour(x_lvl, y_lvl, z_2d
                , c = cgrad(:tempo)
                , clim = (0.35,1), lw = 1.5
                , ylims = (0,41), xlims = (0,21)
                , leg = false
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                )
            scatter!(tmp[!,:P], tmp[!,:S] .- tmp[!,:P], mz = tmp[!,:Pin]
                ,c = cgrad(:tempo), clim = (0.35,1), msw = 0, ma = 0.3, ms = 3
                )
            annotate!(2, 38, text(lti, fntlg))
            d == 1 ? xaxis!("Producer richness") : Nothing
            shape == "Cascade" ? yaxis!("Consumer richness") : Nothing
            #annotate!(2,38, text(lti, 13))
            xticks!([0:5:20;], string.([0:5:20;]))
            yticks!([0:10:40;], string.([0:10:40;]))
            push!(p, thisp)

            if d == 3
                if shape == "Cascade"
                    title!("CASCADE", titlefont = fntlg)
                else
                    title!("PYRAMID", titlefont = fntlg)
                end
            end

            b = DataFrame(P = tmp[!,:bP], H = tmp[!,:bH], C = tmp[!,:bO] .+ tmp[!,:bC])
            pj = plot(legend=false,grid=false,foreground_color_border=:white)
            lmar = shape == "Cascade" ? 15Plots.mm : 7Plots.mm
            for c in 1:size(b,2)
                bm = mean(b[:,c])
                plgn = Shape([-bm/2,bm/2,bm/2,-bm/2], [0+c,0+c,1+c,1+c])
                plot!(plgn, fill = (0, 0.3, :grey), linecolor = :grey, left_margin = lmar, xlims = (-4,4), ylims = (1,5))
            end
            annotate!(0, 4.5, text("N = $N", fntlg))
            xaxis!(false)
            yticks!([1.5:1:3.5;], ["P", "C1", "C2"])
            if shape == "Cascade"
                if d == 1
                    annotate!(-8, 2.5, text("BOTTOM HEAVY", fntlg_r))
                elseif d == 2
                    annotate!(-8, 2.5, text("MIDDLE HEAVY", fntlg_r))
                else
                    annotate!(-8, 2.5, text("TOP HEAVY", fntlg_r))
                end
            end
            push!(p, pj)

        else

            #legend
            pj = plot([NaN],[NaN],grid=false
                , mz = fw[!,:Z], clim = (0.35,1)
                , c = cgrad(:tempo)
                , xaxis = false, yaxis = false
                , xlims = (0,1), ylims = (0,1)
                , legend = :best, label = ""
                #, legendfont = fntlg
                , legendfontsize = 12
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                , colorbar_title = "Prod. intake"
                , foreground_color_legend = nothing)
            push!(p, pj)
            pempty = plot([NaN],[NaN],grid=false,leg=false,axis=false)
            push!(p, pempty)
        end
    end
end
l = grid(3,4, widths = [.15,.35,.35,.15])
porder = p[vcat(2,1,3,4,6,5,7,8,10,9,11,12)]
plot(porder..., layout = l, size = (1000,970), bottom_margin = 15Plots.mm)
savefig("appendices/figures/contour_productivity.png")

#FIGURE S - contour herbivory
p = []
ltr1 = ["A", "B", "C", "D", "E"]
global i = 0
for d in [3,2,1]
    for shape in ["Cascade", "Pyramid"]
        toplot = (fw[!,:shape] .== [shape]) .& (dom .== d)
        #toplot = (dom .== d)
        N = sum(toplot)
        if sum(toplot) != 0
            #get letters for label
            global i = i .+ 1
            lti = ltr1[i]
            #data subset
            tmp = fw[vec(toplot),:]

            tmpdf = DataFrame(P = round.(tmp[!,:P] ./ 10) .* 10
                , A = round.((tmp[!,:S] .- tmp[!,:P]) ./ 5) .* 5
                , B = tmp[!,:PHin] .+ tmp[!,:POin])
            b_mean = aggregate(tmpdf, [:A, :P], mean, skipmissing = true)

            x = b_mean[!,:P]
            y = b_mean[!,:A]
            z = b_mean[!,:B_mean]
            x_lvl = sort(unique(x))
            y_lvl = sort(unique(y))
            z_2d = zeros(length(y_lvl), length(x_lvl))
            for m in 1:length(y_lvl)
                for l in 1:length(x_lvl)
                    xval = x_lvl[l]:x_lvl[l] + 5
                    yval = y_lvl[m]:y_lvl[m] + 10
                    tmpz = [k ∈ xval for k in x] .& [k ∈ yval for k in y]
                    z_val = z[tmpz]
                    z_2d[m,l] = mean(z_val)
                end
            end
            thisp = contour(x_lvl, y_lvl, z_2d
                , c = cgrad(:greys_r)
                , clim = (0,0.4), lw = 1.5
                , ylims = (0,41), xlims = (0,21)
                , leg = false
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                )
            scatter!(tmp[!,:P], tmp[!,:S] .- tmp[!,:P], mz = tmp[!,:PHin] .+ tmp[!,:POin]
                ,c = cgrad(:greys_r), clim = (0,0.4), msw = 0, ma = 0.3, ms = 3
                )
            annotate!(2, 38, text(lti, fntlg))
            d == 1 ? xaxis!("Producer richness") : Nothing
            shape == "Cascade" ? yaxis!("Consumer richness") : Nothing
            #annotate!(2,38, text(lti, 13))
            xticks!([0:5:20;], string.([0:5:20;]))
            yticks!([0:10:40;], string.([0:10:40;]))
            push!(p, thisp)

            if d == 3
                if shape == "Cascade"
                    title!("CASCADE", titlefont = fntlg)
                else
                    title!("PYRAMID", titlefont = fntlg)
                end
            end

            b = DataFrame(P = tmp[!,:bP], H = tmp[!,:bH], C = tmp[!,:bO] .+ tmp[!,:bC])
            pj = plot(legend=false,grid=false,foreground_color_border=:white)
            lmar = shape == "Cascade" ? 15Plots.mm : 7Plots.mm
            for c in 1:size(b,2)
                bm = mean(b[:,c])
                plgn = Shape([-bm/2,bm/2,bm/2,-bm/2], [0+c,0+c,1+c,1+c])
                plot!(plgn, fill = (0, 0.3, :grey), linecolor = :grey, left_margin = lmar, xlims = (-4,4), ylims = (1,5))
            end
            annotate!(0, 4.5, text("N = $N", fntlg))
            xaxis!(false)
            yticks!([1.5:1:3.5;], ["P", "C1", "C2"])
            if shape == "Cascade"
                if d == 1
                    annotate!(-8, 2.5, text("BOTTOM HEAVY", fntlg_r))
                elseif d == 2
                    annotate!(-8, 2.5, text("MIDDLE HEAVY", fntlg_r))
                else
                    annotate!(-8, 2.5, text("TOP HEAVY", fntlg_r))
                end
            end
            push!(p, pj)

        else

            #legend
            pj = plot([NaN],[NaN],grid=false
                , mz = fw[!,:Z], clim = (0,0.4)
                , c = cgrad(:greys_r)
                , xaxis = false, yaxis = false
                , xlims = (0,1), ylims = (0,1)
                , legend = :best, label = ""
                #, legendfont = fntlg
                , legendfontsize = 12
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                , colorbar_title = "Herbivory"
                , foreground_color_legend = nothing)
            push!(p, pj)
            pempty = plot([NaN],[NaN],grid=false,leg=false,axis=false)
            push!(p, pempty)
        end
    end
end
l = grid(3,4, widths = [.15,.35,.35,.15])
porder = p[vcat(2,1,3,4,6,5,7,8,10,9,11,12)]
plot(porder..., layout = l, size = (1000,970), bottom_margin = 15Plots.mm)
savefig("appendices/figures/contour_herbivory.png")

#FIGURE S - contour secondary consumption
p = []
ltr1 = ["A", "B", "C", "D", "E"]
global i = 0
for d in [3,2,1]
    for shape in ["Cascade", "Pyramid"]
        toplot = (fw[!,:shape] .== [shape]) .& (dom .== d)
        #toplot = (dom .== d)
        N = sum(toplot)
        if sum(toplot) != 0
            #get letters for label
            global i = i .+ 1
            lti = ltr1[i]
            #data subset
            tmp = fw[vec(toplot),:]

            tmpdf = DataFrame(P = round.(tmp[!,:P] ./ 10) .* 10
                , A = round.((tmp[!,:S] .- tmp[!,:P]) ./ 5) .* 5
                , B = tmp[!,:OOin] .+ tmp[!,:OCin] .+ tmp[!,:CCin] .+ tmp[!,:COin] .+ tmp[!,:HCin])
            b_mean = aggregate(tmpdf, [:A, :P], mean, skipmissing = true)

            x = b_mean[!,:P]
            y = b_mean[!,:A]
            z = b_mean[!,:B_mean]
            x_lvl = sort(unique(x))
            y_lvl = sort(unique(y))
            z_2d = zeros(length(y_lvl), length(x_lvl))
            for m in 1:length(y_lvl)
                for l in 1:length(x_lvl)
                    xval = x_lvl[l]:x_lvl[l] + 5
                    yval = y_lvl[m]:y_lvl[m] + 10
                    tmpz = [k ∈ xval for k in x] .& [k ∈ yval for k in y]
                    z_val = z[tmpz]
                    z_2d[m,l] = mean(z_val)
                end
            end
            thisp = contour(x_lvl, y_lvl, z_2d
                , c = cgrad(:plasma_r)
                , clim = (0,0.4), lw = 1.5
                , ylims = (0,41), xlims = (0,21)
                , leg = false
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                )
            scatter!(tmp[!,:P], tmp[!,:S] .- tmp[!,:P], mz = tmp[!,:PHin] .+ tmp[!,:POin]
                ,c = cgrad(:plasma_r), clim = (0,0.4), msw = 0, ma = 0.3, ms = 3
                )
            annotate!(2, 38, text(lti, fntlg))
            d == 1 ? xaxis!("Producer richness") : Nothing
            shape == "Cascade" ? yaxis!("Consumer richness") : Nothing
            #annotate!(2,38, text(lti, 13))
            xticks!([0:5:20;], string.([0:5:20;]))
            yticks!([0:10:40;], string.([0:10:40;]))
            push!(p, thisp)

            if d == 3
                if shape == "Cascade"
                    title!("CASCADE", titlefont = fntlg)
                else
                    title!("PYRAMID", titlefont = fntlg)
                end
            end

            b = DataFrame(P = tmp[!,:bP], H = tmp[!,:bH], C = tmp[!,:bO] .+ tmp[!,:bC])
            pj = plot(legend=false,grid=false,foreground_color_border=:white)
            lmar = shape == "Cascade" ? 15Plots.mm : 7Plots.mm
            for c in 1:size(b,2)
                bm = mean(b[:,c])
                plgn = Shape([-bm/2,bm/2,bm/2,-bm/2], [0+c,0+c,1+c,1+c])
                plot!(plgn, fill = (0, 0.3, :grey), linecolor = :grey, left_margin = lmar, xlims = (-4,4), ylims = (1,5))
            end
            annotate!(0, 4.5, text("N = $N", fntlg))
            xaxis!(false)
            yticks!([1.5:1:3.5;], ["P", "C1", "C2"])
            if shape == "Cascade"
                if d == 1
                    annotate!(-8, 2.5, text("BOTTOM HEAVY", fntlg_r))
                elseif d == 2
                    annotate!(-8, 2.5, text("MIDDLE HEAVY", fntlg_r))
                else
                    annotate!(-8, 2.5, text("TOP HEAVY", fntlg_r))
                end
            end
            push!(p, pj)

        else

            #legend
            pj = plot([NaN],[NaN],grid=false
                , mz = fw[!,:Z], clim = (0,0.4)
                , c = cgrad(:plasma_r)
                , xaxis = false, yaxis = false
                , xlims = (0,1), ylims = (0,1)
                , legend = :best, label = ""
                #, legendfont = fntlg
                , legendfontsize = 12
                , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg
                , colorbar_title = "Sec. consumption"
                , foreground_color_legend = nothing)
            push!(p, pj)
            pempty = plot([NaN],[NaN],grid=false,leg=false,axis=false)
            push!(p, pempty)
        end
    end
end
l = grid(3,4, widths = [.15,.35,.35,.15])
porder = p[vcat(2,1,3,4,6,5,7,8,10,9,11,12)]
plot(porder..., layout = l, size = (1000,970), bottom_margin = 15Plots.mm)
savefig("appendices/figures/contour_consumption.png")

#FIGURE 2 - characteristics of different trophic structure
tmpdf = DataFrame(d = vec(dom), s = fw[!,:shape], h = round.(fw[!,:Height_mean], digits = 1), C = fw[!,:L] ./ (fw[!,:S] .^2), Z = log10.(fw[!,:Z]))
m = DataFrame(d = vec(dom), s = fw[!,:shape], IGP = fw[!,:OMN] ./ fw[!,:TOT], LFC = fw[!,:LFC] ./ fw[!,:TOT], AC = fw[!,:AC] ./ fw[!,:TOT], EC = fw[!,:EC] ./ fw[!,:TOT], C = fw[!,:L] ./ (fw[!,:S] .^2))

cascades = tmpdf[tmpdf[!,:s] .== "Cascade",:]
pyramids = tmpdf[tmpdf[!,:s] .== "Pyramid",:]

pz = @df cascades violin(:d,:Z, side=:right, fc = RGB(0.45, 0.45, 0.45), leg = false
    , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg, left_margin = 5Plots.mm)
xticks!([1:1:3;], ["BH", "MH", "TH"])
@df pyramids violin!(:d,:Z, side=:left, fc = RGB(0.75, 0.75, 0.75))
df = aggregate(tmpdf[!,[:d, :Z]], :d, mean)
@df df scatter!(:d, :Z_mean, ms = 5, mc = :black, msc = :white)
yaxis!("Average Z")
annotate!(0.75, 4.5, text("B", 14))
#savefig("manuscript/figures/Z.png")

ph = @df cascades violin(:d,:h, side=:right, fc = RGB(0.45, 0.45, 0.45), leg = false
    , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg, left_margin = 5Plots.mm)
xticks!([1:1:3;], ["BH", "MH", "TH"])
@df pyramids violin!(:d,:h, side=:left, fc = RGB(0.75, 0.75, 0.75))
df = aggregate(tmpdf[!,[:d, :h]], :d, mean)
@df df scatter!(:d, :h_mean, ms = 5, mc = :black, msc = :white)
yaxis!("Average height")
annotate!(0.75, 4.1, text("A", 14))
#savefig("manuscript/figures/height.png")

pc = @df cascades violin(:d,:C, side=:right, fc = RGB(0.45, 0.45, 0.45), leg = false
    , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg, left_margin = 5Plots.mm)
xticks!([1:1:3;], ["BH", "MH", "TH"])
@df pyramids violin!(:d,:C, side=:left, fc = RGB(0.75, 0.75, 0.75),)
df = aggregate(tmpdf[!,[:d, :C]], :d, mean)
@df df scatter!(:d, :C_mean, ms = 5, mc = :black, msc = :white)
yaxis!("Connectance")
annotate!(0.75, 0.35, text("C", 14))
#savefig("manuscript/figures/connectance.png")

m_cascades = m[m[!,:s] .== "Cascade",:]
m_pyramids = m[m[!,:s] .== "Pyramid",:]

p1 = @df m_cascades violin(:d,:IGP, side=:right, fc = RGB(0.45, 0.45, 0.45), leg = false, ylims = (0,1)
    , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg, left_margin = 5Plots.mm)
xticks!([1:1:3;], ["BH", "MH", "TH"])
@df m_pyramids violin!(:d,:IGP, side=:left, fc = RGB(0.75, 0.75, 0.75))
df = aggregate(m[!,[:d, :IGP]], :d, mean)
@df df scatter!(:d, :IGP_mean, ms = 5, mc = :black, msc = :white)
yaxis!("Freq. omnivory / IGP")
annotate!(0.75, 0.95, text("D", 14))

p2 = @df m_cascades violin(:d,:LFC, side=:right, fc = RGB(0.45, 0.45, 0.45), leg = false, ylims = (0,1)
    , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg, left_margin = 5Plots.mm)
xticks!([1:1:3;], ["BH", "MH", "TH"])
@df m_pyramids violin!(:d,:LFC, side=:left, fc = RGB(0.75, 0.75, 0.75))
df = aggregate(m[!,[:d, :LFC]], :d, mean)
@df df scatter!(:d, :LFC_mean, ms = 5, mc = :black, msc = :white)
yaxis!("Freq. food chain")
annotate!(0.75, 0.95, text("E", 14))

p3 = @df m_cascades violin(:d,:AC, side=:right, fc = RGB(0.45, 0.45, 0.45), leg = false, ylims = (0,1)
    , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg, left_margin = 5Plots.mm)
xticks!([1:1:3;], ["BH", "MH", "TH"])
@df m_pyramids violin!(:d,:AC, side=:left, fc = RGB(0.75, 0.75, 0.75))
df = aggregate(m[!,[:d, :AC]], :d, mean)
@df df scatter!(:d, :AC_mean, ms = 5, mc = :black, msc = :white)
yaxis!("Freq. apparent compet.")
annotate!(0.75, 0.95, text("F", 14))

p4 = @df m_cascades violin(:d,:EC, side=:right, fc = RGB(0.45, 0.45, 0.45), leg = false, ylims = (0,1)
    , ytickfont=fntlg, xtickfont = fntlg, guidefont = fntlg, left_margin = 5Plots.mm)
xticks!([1:1:3;], ["BH", "MH", "TH"])
@df m_pyramids violin!(:d,:EC, side=:left, fc = RGB(0.75, 0.75, 0.75),)
df = aggregate(m[!,[:d, :EC]], :d, mean)
yaxis!("Freq. exploitative compet.")
@df df scatter!(:d, :EC_mean, ms = 5, mc = :black, msc = :white)
annotate!(0.75, 0.95, text("G", 14))

#plot(p1,p2,p3,p4, layout = grid(2,2), size = (900,600))
#savefig("manuscript/figures/returnofthemotifs.png")

pe = scatter([NaN], [NaN], xaxis = false, yaxis = false, grid = false, leg = false)
pl = scatter([NaN NaN NaN], [NaN NaN NaN], markershape = [:rect :rect :circle]
    , mc = [RGB(0.45, 0.45, 0.45) RGB(0.75, 0.75, 0.75) :black]
    , labels = ["Cascade-shaped" "Pyramid-shaped" "Average"]
    , msc = :black, xaxis = false, yaxis = false, grid = false
    , legendfont = fntlg
    , legend = :best)

plot(ph, pz, pc, p1, p2, p3, p4, pe, pl, size = (800,800))
savefig("manuscript/figures/structure.png")
