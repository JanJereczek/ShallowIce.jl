include("intro.jl")
using CairoMakie

function main(; N = 101, tend = 2e6)
    data = jldopen(datadir("hysteresis_N=$(N)_tend=$(tend).jld2"))
    ht, bt, t_out, a_out = data["ht"], data["bt"], data["t_out"], data["a_out"]

    nx, nt = size(ht)
    a_f = [a_out[j][1] for j in eachindex(a_out)]
    L_raw = [findlast(ht[:, j] .> 1e-3) for j in axes(ht, 2)]
    L = Int.(L_raw[1:nt-1])

    ctickvec = 0:5e5:2e6
    cticklabel = Int.(ctickvec ./ 1e5)
    ft = 32
    fig = Figure(resolution = (1200, 900), fontsize = ft)
    ax = Axis(
        fig[1, 1],
        xlabel = L"Mass balance at top of ice sheet (m/yr) $\,$",
        ylabel = L"Last ice-covered point (1) $\,$",
    )
    lh = lines!(ax, a_f[1:nt-1], L, color = t_out[1:nt-1], colormap = :jet)
    cticks = (ctickvec, [L"%$(c) $\,$" for c in cticklabel])
    Colorbar(fig[1, 2], lh, label = L"Time elapsed ($10^5$ yr)", ticks = cticks, height = Relative(0.5))
    save(plotsdir("hysteresis.png"), fig)
    save(plotsdir("hysteresis.pdf"), fig)
end
main(N=51)
# SH co2 is leading temperature signal, NH co2 is lagging compared to temperature.