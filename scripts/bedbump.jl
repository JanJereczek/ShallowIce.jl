include("intro.jl")
using CairoMakie

function main(; N = 101, tmax = 60e3)
    N = 101
    L = 1.5e6
    tspan = (0.0, tmax)
    dt = 2.0
    bc = "zero_flow"

    omega = ComputationDomain(L, N, tspan, dt, bc)
    bedhead = 2e3               # m
    bedslope = -0.002           # (m/m)
    bedbump_center = 1e6        # m
    bedbump_height = 500.0      # m
    bedbump_width = 100e3       # m
    h0 = fill(0.0, N)
    b0 = bumpy_bed(omega.xH, bedhead, bedslope, bedbump_center, bedbump_height, bedbump_width)
    
    p = Params(isostasy_on = true, accumulation = linear_accumulation)
    iss = IcesheetState(N, h0, b0)
    dt_out = 100
    t_out = Int.(collect(tspan[1]:dt_out:tspan[2]))
    y_out = fill(1, length(t_out))
    sstruct = SuperStruct(p, omega, iss)
    ht, bt = forward_sia(sstruct, dt_out = dt_out)

    idx = Observable(1)
    stride = 1
    hplot = @lift(ht[:, $idx])
    bplot = @lift(bt[:, $idx])
    splot = @lift((bt .+ ht)[:, $idx])
    yplot = @lift(y_out[1:$idx])
    tplot = @lift("t = $(t_out[$idx]) yr")
    x = Int.(omega.xH .÷ 1000)

    lw = 5
    ft = 32
    fig = Figure(resolution = (1600, 900), fontsize = 32)
    ax = Axis(
        fig[1, 1],
        xlabel = L"Position (km) $\,$",
        ylabel = L"Surface elevation (m) $\,$",
    )
    lh = lines!(ax, x, splot, color = :cornflowerblue,
        linewidth = lw, label = L"z_{ice} $\,$")
    lb = lines!(ax, x, bplot, color = :gray10, linewidth = lw, label = L"b $\,$")
    lb0 = lines!(ax, x, b0, color = :gray50, linestyle = :dash,
        linewidth = lw, label = L"b_{0} $\,$")
    text!(ax, 1300, 3500, text = tplot, fontsize = ft)

    ylims!(ax, (-1500, 4000))
    record(fig, "plots/bumpy_isostasy_N=$(N)_$bc.gif", axes(ht, 2)[1:stride:end]) do i
        idx[] = i
    end
end

main()