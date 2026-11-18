# Probe: the energy that the collision parts of the Fokker-Planck term give to a Maxwellian at the
# background's own (n, T) on the grid of the case, against T; and to a Maxwellian plus a cold
# population in the lowest cell (the shape of the particle source).
include(joinpath(@__DIR__, "common.jl"))
using Serialization
SCHEME = length(ARGS) >= 1 ? Symbol(ARGS[1]) : :tcdf
LOG = length(ARGS) >= 2 ? Symbol(ARGS[2]) : :energy_dependent
if LOG !== :energy_dependent
    @eval collision_settings(compat; mode = :full) = CollisionSettings(; mode, screening = :partial, coulomb_logarithm = $(QuoteNode(LOG)))
end
println("scheme = $SCHEME, coulomb_logarithm = $LOG")
case = build(START; compat = :runes2p, runaway = false, scheme = SCHEME); pc = proto_case(case)
starts = deserialize(joinpath(OUT, "starts.jls"))
u = copy(starts.ub); t = TIMES[2]
fp = fp_of(pc)
n = vals(pc, u, :cold_electron_density)
println("rates [1/s] = power / (3/2 e n T); +: f gains energy. Columns r1 r2 (ions of the start state at both radii)")
println("  T [eV] | e-e (−L_fe) on M | sinks (−L_s) on M | e-e on M + 1% of n in the lowest cell | sinks on that")
for T in (1500.0, 500.0, 139.0, 50.0, 10.0)
    s = setvals!(pc, copy(u), :cold_electron_temperature, fill(T, pc.nr))
    F = maxwellian_array(pc, n, fill(T, pc.nr))
    s = withf(pc, s, F)
    prepare_terms!(pc, s, t)
    W = THERMAL_ENERGY_FACTOR .* n .* T
    fe = -tv(pc, KineticEnergyLossTerm(fp; part = :free_electrons, coefficient = 1.0), s, t) ./ W
    sk = -tv(pc, KineticEnergyLossTerm(fp; part = :sinks, coefficient = 1.0), s, t) ./ W
    F2 = copy(F)
    for ir in 1:pc.nr
        D, _ = weights(pc.grid, ir)
        F2[1, :, ir] .+= 0.01 * n[ir] / sum(D[1, :])
    end
    s2 = withf(pc, s, F2)
    fe2 = -tv(pc, KineticEnergyLossTerm(fp; part = :free_electrons, coefficient = 1.0), s2, t) ./ W
    sk2 = -tv(pc, KineticEnergyLossTerm(fp; part = :sinks, coefficient = 1.0), s2, t) ./ W
    println(@sprintf("  %6.0f | %s | %s | %s | %s", T, v(fe), v(sk), v(fe2), v(sk2)))
end
println("END")
