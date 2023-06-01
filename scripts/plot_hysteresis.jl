include("intro.jl")
using CairoMakie

function main(; N = 101, tend = 2e6, isostasy_on = false)
    data = jldopen(datadir("hysteresis_N=$(N)_tend=$(tend)_isostasy=$(isostasy_on).jld2"))
    ht, bt, t_out, z_out = data["ht"], data["bt"], data["t_out"], data["z_out"]

    nx, nt = size(ht)
    L_raw = [findlast(ht[:, j] .> 1e-3) for j in axes(ht, 2)]
    L = copy(L_raw)
    L[L_raw .== nothing] .= 0
    L = Int.(L)

    ctickvec = 0:5e5:2e6
    cticklabel = Int.(ctickvec ./ 1e5)
    cticks = (ctickvec, [L"%$(c) $\,$" for c in cticklabel])
    ft = 32
    fig = Figure(resolution = (1600, 900), fontsize = ft)
    
    axf = Axis(
        fig[1, 1],
        ylabel = L"Equilibrium line (m) $\,$",
        xticklabelsvisible = false,
    )
    axl = Axis(
        fig[2, 1],
        xlabel = L"Time ($10^5$ yr)",
        ylabel = L"Equilibrium line (m) $\,$",
        xticks = cticks,
    )
    axh = Axis(
        fig[:, 2],
        xlabel = L"Equilibrium line (m) $\,$",
        ylabel = L"Last ice-covered point (1) $\,$",
    )

    lf = lines!(axf, t_out, z_out)
    ll = lines!(axl, t_out[1:nt], L)
    lh = lines!(axh, z_out[1:nt], L, color = t_out[1:nt], colormap = :jet)
    Colorbar(fig[:, 3], lh, label = L"Time elapsed ($10^5$ yr)", ticks = cticks, height = Relative(0.5))
    save(plotsdir("hysteresis_isostasy=$(isostasy_on).png"), fig)
    save(plotsdir("hysteresis_isostasy=$(isostasy_on).pdf"), fig)
end
main(N=76)