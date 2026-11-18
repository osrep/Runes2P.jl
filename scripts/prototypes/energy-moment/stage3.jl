# Stage 3: (1) the cold row from start (b); (2) the integral of the e-e heating; (3) the cooling alone on a Maxwellian.
include(joinpath(@__DIR__, "common.jl"))
using Serialization
starts = deserialize(joinpath(OUT, "starts.jls"))
base_r = base_case(:runes2p); pc = proto_case(base_r)

# (1) cold row from start (b), state copied by unknown name
u = zeros(length(base_r.layout))
for unk in base_r.layout.unknowns
    name = unk.name
    src = reduce(vcat, collect.(Runes2P._unknown_ranges(pc.layout, unknown_id(pc.layout, name))))
    dst = reduce(vcat, collect.(Runes2P._unknown_ranges(base_r.layout, unknown_id(base_r.layout, name))))
    u[dst] .= starts.ub[src]
end
times = collect(range(0.0, 30e-7; length = 31))
system = EquationSystem(base_r.layout, base_r.equations; time_scale = 1e-7, initial_state = u)
sol = integrate(system, u, times; stepper = FixedSteps())
println("cold row from start (b): ", stats_line(sol))
r30 = deserialize(joinpath(OUT, "results_tcdf_30_pa-pb-cold.jls"))
_, recs_pb, _ = r30[:pb]; _, recs_c, _ = r30[:cold]
println("  t | T cold row from (b) | T (ii) | T (iii) from DREAM's state | n_hot cold row from (b) | n_hot (ii) | n_hot (iii) | E cold(b) | E (ii)")
for k in (1, 2, 4, 11, 21, 31)
    s = sol.states[:, k]
    println(@sprintf("  %.1e | %s | %s | %s | %s | %s | %s | %s | %s", times[k], v(vals(base_r, s, :cold_electron_temperature)), v(recs_pb[k].T), v(recs_c[k].T),
        v(vals(base_r, s, :hot_electron_density)), v(recs_pb[k].n_hot), v(recs_c[k].n_hot), v(vals(base_r, s, :electric_field)), v(recs_pb[k].E)))
end
Tcb = vals(base_r, sol.states[:, 11], :cold_electron_temperature)
pr("(T(ii) − T_coldrow(b))/T_coldrow(b) at 1e-6 s", (recs_pb[11].T .- Tcb) ./ Tcb)
Tcb30 = vals(base_r, sol.states[:, 31], :cold_electron_temperature)
pr("(T(ii) − T_coldrow(b))/T_coldrow(b) at 3e-6 s", (recs_pb[31].T .- Tcb30) ./ Tcb30)

# (2) the integral of the e-e heating and of the sinks over the first 10 steps, against the excess energy of the bulk
for (file, key, label) in (("budgets_tcdf_10_pa-pb-pas-cold-total-dream.jls", :pb, "tcdf (ii)"), ("budgets_exponential_fitting_30_pb-cold.jls", :pb, "exponential fitting (ii)"))
    b = deserialize(joinpath(OUT, file))[key]
    ee = sum(r.ee for r in b[1:10]) .* 1e-7; sk = sum(r.sinks for r in b[1:10]) .* 1e-7
    eeM = sum(r.eeM for r in b[1:10]) .* 1e-7; cool = sum(r.cool for r in b[1:10]) .* 1e-7
    println("  $label, first 10 steps:")
    pr("    ∫ e-e heating dt [J/m^3]", ee); pr("    ∫ sinks dt [J/m^3]", sk); pr("    ∫ e-e on a Maxwellian at (n,T) dt [J/m^3]", eeM)
    pr("    ∫ cooling (P_fluid) dt [J/m^3]", cool)
end
pr("W_below(ii, tcdf) − W_cold(iii) at 1e-6 s [J/m^3]", recs_pb[11].W_below .- recs_c[11].W_cold)
ef = deserialize(joinpath(OUT, "results_exponential_fitting_30_pb-cold.jls"))
pr("W_below(ii, exp. fitting) − W_cold(iii, exp. fitting) at 1e-6 s [J/m^3]", ef[:pb][2][11].W_below .- ef[:cold][2][11].W_cold)

# (3) the cooling alone on a Maxwellian at r1: 30 explicit steps of 1e-7 s at a = 1.15e5 1/s
g = base_r.grid; region = region_of(base_r)
for (T0, a) in ((1580.0, 1.15e5), (139.0, 5e5))
    ir = 1
    F = maxwellian_array(base_r, fill(6.7e19, base_r.nr), fill(T0, base_r.nr))[:, :, ir]
    D, E = weights(g, ir)
    Tm = sum(E .* F) / (THERMAL_ENERGY_FACTOR * sum(D .* F))
    for _ in 1:30
        w = envelope(region, g, ir, Tm)
        df, _ = unit_advection(g, ir, F, w, -1.0)
        F .+= 1e-7 * a .* df
        Tm = sum((1 .- w) .* E .* F) / (THERMAL_ENERGY_FACTOR * sum((1 .- w) .* D .* F))
    end
    n = sum(D .* F)
    FM = maxwellian_array(base_r, fill(n, base_r.nr), fill(Tm, base_r.nr))[:, :, ir]
    # the Maxwellian at the moment temperature has a different moment temperature on the grid; compare also with the one of equal moment
    L1 = sum(D .* abs.(F .- FM)) / n
    println(@sprintf("  cooling alone, T0=%.0f eV, a=%.2e 1/s, 30 steps: T_moment=%.4g eV (%.3f of T0), min f=%.3e, L1 to Maxwellian(n, T_moment)=%.3e",
        T0, a, Tm, Tm / T0, minimum(F), L1))
    m = momentum_grid(g, ir)
    println("   f/f_M: ", join([@sprintf("%.3g:%.3g", m.p_centers[i], F[i, 1] / FM[i, 1]) for i in 1:momentum_cell_count(m) if m.p_centers[i] < 4 * sqrt(2Tm / MC2)], " "))
end
println("END")
