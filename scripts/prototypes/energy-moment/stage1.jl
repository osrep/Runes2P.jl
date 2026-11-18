# Stage 1: static checks on the start states (a) and (b); point 3 at the start; the grid.
include(joinpath(@__DIR__, "common.jl"))

case = base_case(:runes2p)
pc = proto_case(case)
g = case.grid
m1 = momentum_grid(g, 1)
println("p_faces = ", round.(m1.p_faces; sigdigits = 4))
println("xi_faces = ", round.(m1.xi_faces; sigdigits = 4))
println("TIMES = ", TIMES)
ε1 = Runes2P._kinetic_energy(m1.p_centers[1]) / e
println(@sprintf("first cell centre p=%.4g, ε=%.4g eV, floor of T from a bulk in cell 1 = 2/3 ε = %.4g eV", m1.p_centers[1], ε1, 2ε1 / 3))

# quadrature of a Maxwellian on the grid against T
println("\nMaxwellian on the grid: n_grid/n, W_grid/(3/2 n T), T_moment/T (whole grid)")
for T in (2000.0, 1000.0, 500.0, 300.0, 200.0, 139.0, 100.0, 70.0, 50.0, 30.0, 20.0, 10.0, 5.0)
    F = maxwellian_array(case, fill(1e19, case.nr), fill(T, case.nr))
    D, E = weights(g, 1)
    n = sum(D .* F[:, :, 1]); W = sum(E .* F[:, :, 1])
    println(@sprintf("  T=%7.1f eV  n_grid/n=%.4f  W/(1.5nT)=%.4f  T_mom/T=%.4f", T, n / 1e19, W / (1.5e * 1e19 * T), W / (1.5e * n * T)))
end

# the below + above = whole check
u0 = initial_state(pc, START)
solve_rows!(pc, u0, DERIVED)
t1 = TIMES[2]
row = radial_row(pc)
Wb = tv(pc, EnergyMomentTerm(:f, g; region = pc.below), u0, t1)
Wh = tv(pc, EnergyMomentTerm(:f, g; region = region_of(pc)), u0, t1)
Wa = tv(pc, EnergyMomentTerm(:f, g), u0, t1)
pr("below + hot − whole [J/m^3]", Wb .+ Wh .- Wa)
TD = vals(pc, u0, :cold_electron_temperature)
pr("DREAM T_cold [eV]", TD)

# start (a): scan of F(T)/(1.5 e n T) and the root
println("\nstart (a): scan of the moment row F(T) = M_below(T) − 3/2 e n_cold(T) T, relative to 3/2 e n T")
for x in (0.5, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 4.0, 6.0)
    F, s = Frow(pc, u0, x .* TD)
    n = vals(pc, s, :cold_electron_density)
    println(@sprintf("  T/T_D=%.2f  F/(1.5enT) = [%s]   p0 = [%s]", x, v(F ./ (THERMAL_ENERGY_FACTOR .* n .* x .* TD)),
        v(p0_of.(Ref(region_of(pc)), x .* TD))))
end
Ta, ua = moment_temperature(pc, u0, 0.9 .* TD, 6.0 .* TD)
pr("start (a) T from the moment [eV]", Ta); pr("T_a / T_D", Ta ./ TD)
# amplification: (dM_below/dT)/(3/2 e n) at the root, at fixed f
h = 1e-4 .* Ta
Mp = moments(pc, setvals!(pc, copy(ua), :cold_electron_temperature, Ta .+ h)).Wbelow
Mm = moments(pc, setvals!(pc, copy(ua), :cold_electron_temperature, Ta .- h)).Wbelow
na = vals(pc, ua, :cold_electron_density)
pr("(dM_below/dT)/(3/2 e n) at the root [-]", (Mp .- Mm) ./ (2h) ./ (THERMAL_ENERGY_FACTOR .* na))
pr("n_cold(a) − density below [m^-3]", na .- moments(pc, ua).Nbelow)

# start (b): bulk replaced by c M(T*) (1 + α ξ p) below the threshold at T_D, density, moment-T and bulk current kept
function start_b(pc, u0, TD)
    F0 = farr(pc, u0); M0 = fmoments(pc, F0, TD)
    u = copy(u0); solve_rows!(pc, u, DERIVED)
    jb0 = tv(pc, CurrentMomentTerm(:f, g; pitch_rule = :cell_center), u, t1) .- vals(pc, u, :hot_electron_current_density)
    Tstar = copy(TD); F = copy(F0); info = Dict{Symbol, Any}()
    for ir in 1:pc.nr
        m = momentum_grid(g, ir); D, E = weights(g, ir)
        w = envelope(region_of(pc), g, ir, TD[ir])
        ξ = m.xi_centers; p = m.p_centers
        mk(M, c) = [(1 - w[i]) * c * M[i, j] + w[i] * F0[i, j, ir] for i in eachindex(p), j in eachindex(ξ)]
        function fb(Ts)
            M = [exp(-Runes2P._kinetic_energy(p[i]) / (Runes2P.Constants.ELEMENTARY_CHARGE * Ts)) for i in eachindex(p), j in eachindex(ξ)]
            # c: keep the density of f
            base = sum(D .* mk(M, 0.0)); unit = sum(D .* (1 .- w) .* M)
            c = (M0.Nall[ir] - base) / unit
            return mk(M, c), c
        end
        # G(T*) = M_below(f_b; T_D) − 3/2 e n_cold(f_b) T_D with n_cold = n_free − n_hot(f_b) − n_re
        nfree = vals(pc, u, :free_electron_density)[ir]; nre = vals(pc, u, :runaway_density)[ir]
        G(Ts) = (f = fb(Ts)[1]; sum((1 .- w) .* E .* f) - THERMAL_ENERGY_FACTOR * (nfree - sum(w .* D .* f) - nre) * TD[ir])
        lo, hi = 0.5TD[ir], 2TD[ir]
        for x in (0.3, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0)
            ff, cc = fb(x * TD[ir])
            println(@sprintf("   r%d T*/T_D=%.2f G/(1.5 e n T)=%.4e c=%.4e Mbelow=%.4e nhot=%.4e", ir, x, G(x * TD[ir]) / (THERMAL_ENERGY_FACTOR * nfree * TD[ir]), cc,
                sum((1 .- w) .* E .* ff), sum(w .* D .* ff)))
        end
        G(lo) * G(hi) < 0 || error("no bracket for T* at radius $ir")
        for _ in 1:80
            mid = (lo + hi) / 2
            (G(mid) * G(lo) > 0) ? (lo = mid) : (hi = mid)
        end
        Tstar[ir] = (lo + hi) / 2
        f, c = fb(Tstar[ir])
        # α: keep the bulk current j(f) − j_hot(f)
        F[:, :, ir] .= f
        info[Symbol("c$ir")] = c
    end
    ub = withf(pc, u, F); solve_rows!(pc, ub, DERIVED)
    jb1 = tv(pc, CurrentMomentTerm(:f, g; pitch_rule = :cell_center), ub, t1) .- vals(pc, ub, :hot_electron_current_density)
    # current of the unit α perturbation (linear in α)
    Fα = copy(F)
    for ir in 1:pc.nr
        m = momentum_grid(g, ir); w = envelope(region_of(pc), g, ir, TD[ir])
        M = [exp(-Runes2P._kinetic_energy(m.p_centers[i]) / (e * Tstar[ir])) for i in eachindex(m.p_centers), j in eachindex(m.xi_centers)]
        c = info[Symbol("c$ir")]
        Fα[:, :, ir] .+= [(1 - w[i]) * c * M[i, j] * m.xi_centers[j] * m.p_centers[i] for i in 1:size(M, 1), j in 1:size(M, 2)]
    end
    uα = withf(pc, ub, Fα); solve_rows!(pc, uα, DERIVED)
    jα = tv(pc, CurrentMomentTerm(:f, g; pitch_rule = :cell_center), uα, t1) .- vals(pc, uα, :hot_electron_current_density) .- jb1
    α = (jb0 .- jb1) ./ jα
    Ffinal = copy(F)
    for ir in 1:pc.nr
        Ffinal[:, :, ir] .+= α[ir] .* (Fα[:, :, ir] .- F[:, :, ir])
    end
    ub = withf(pc, u, Ffinal); solve_rows!(pc, ub, DERIVED)
    return ub, Tstar, α, jb0
end
ub, Tstar, α, jb0 = start_b(pc, u0, TD)
pr("start (b): Maxwellian temperature T* [eV]", Tstar); pr("T*/T_D", Tstar ./ TD)
pr("start (b): α [1/(m_e c)]", α)
Fb = farr(pc, ub); pr("start (b): min f", [minimum(Fb[:, :, ir]) for ir in 1:pc.nr])
Tb_check, _ = moment_temperature(pc, ub, 0.5 .* TD, 2.0 .* TD)
pr("start (b): T from the moment / T_D", Tb_check ./ TD)
pr("start (b): bulk current before [A/m^2]", jb0)
pr("start (b): bulk current after [A/m^2]", tv(pc, CurrentMomentTerm(:f, g; pitch_rule = :cell_center), ub, t1) .- vals(pc, ub, :hot_electron_current_density))
pr("start (b): ∫f − ∫f_DREAM [m^-3]", moments(pc, ub).Nall .- moments(pc, u0).Nall)

# point 3 at the starts: parts of the Fokker-Planck term, whole grid
function fp_parts(pc, u, t, label)
    prepare_terms!(pc, u, t)
    fp = fp_of(pc)
    W, N = term_f(pc, fp, u, t)
    Lfe = tv(pc, KineticEnergyLossTerm(fp; part = :free_electrons, coefficient = 1.0), u, t)
    Ls = tv(pc, KineticEnergyLossTerm(fp; part = :sinks, coefficient = 1.0), u, t)
    Lsb = tv(pc, KineticEnergyLossTerm(fp; part = :sinks, region = pc.below, coefficient = 1.0), u, t)
    Lfeb = tv(pc, KineticEnergyLossTerm(fp; part = :free_electrons, region = pc.below, coefficient = 1.0), u, t)
    out = tv(pc, MomentumBoundaryFluxTerm(:f, fp; weight = :energy), u, t)
    println("  --- $label")
    pr("FP energy moment, whole grid [W/m^3]", W)
    pr("free-electron loss L_fe whole grid [W/m^3] (+: f loses)", Lfe)
    pr("free-electron loss below the threshold [W/m^3]", Lfeb)
    pr("sinks L_s whole grid [W/m^3]", Ls); pr("sinks below the threshold [W/m^3]", Lsb)
    pr("outflow energy at p_max (term value) [W/m^3]", out)
    pr("field part = W + L_fe + L_s + outflow [W/m^3]", W .+ Lfe .+ Ls .+ out)
    pr("E (j(f)) [W/m^3]", vals(pc, u, :electric_field) .* tv(pc, CurrentMomentTerm(:f, g; pitch_rule = :cell_center), u, t))
    # the fluid channels
    c = pc.cooling; c.key = nothing; compute!(c, u, pc.layout, t)
    for k in sort(collect(keys(c.parts)))
        pr("P_fluid channel: $k [W/m^3]", c.parts[k])
    end
    pr("P_fluid total [W/m^3]", c.P); pr("cooling a_k [1/s]", c.a); pr("cooling rate max |a| p g J/(J Δp) [1/s]", c.rate)
    pr("E j_ohm [W/m^3]", vals(pc, u, :electric_field) .* vals(pc, u, :ohmic_current_density))
    # L_fe on a Maxwellian at (n_cold, T) of the state
    T = vals(pc, u, :cold_electron_temperature); n = vals(pc, u, :cold_electron_density)
    uM = withf(pc, u, maxwellian_array(pc, n, T))
    pr("L_fe on a Maxwellian at (n_cold,T) [W/m^3]", tv(pc, KineticEnergyLossTerm(fp; part = :free_electrons, coefficient = 1.0), uM, t))
    pr("  relative to 3/2 e n T [1/s]", tv(pc, KineticEnergyLossTerm(fp; part = :free_electrons, coefficient = 1.0), uM, t) ./ (1.5e .* n .* T))
    pr("L_fe / (3/2 e n T) [1/s]", Lfe ./ (1.5e .* n .* T))
    return nothing
end
println("\npoint 3 at the starts (fp prepared at the start state, t = t1)")
fp_parts(pc, ua, t1, "start (a)")
fp_parts(pc, ub, t1, "start (b)")
u0c = copy(u0); solve_rows!(pc, u0c, DERIVED)
fp_parts(pc, u0c, t1, "DREAM's state as is, T = T_D (not on the moment row)")

using Serialization
serialize(joinpath(OUT, "starts.jls"), (; ua, ub, Ta, Tstar, α))
println("END")
