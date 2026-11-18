# Stage 2: the runs and the comparison. ARGS[1] = number of steps of 1e-7 s (10: DREAM's window).
include(joinpath(@__DIR__, "common.jl"))
using Serialization
const NSTEPS = parse(Int, ARGS[1])
const VARIANTS = length(ARGS) >= 2 ? Symbol.(split(ARGS[2], ",")) : [:pa, :pb, :pas, :cold, :total, :dream]
times = NSTEPS <= 10 ? TIMES[1:(NSTEPS + 1)] : collect(range(0.0, NSTEPS * 1e-7; length = NSTEPS + 1))
println("times: ", NSTEPS, " steps to ", times[end])

base_r = base_case(:runes2p)
pc = proto_case(base_r)
pcs = proto_case(base_r; compensate_sinks = true)
starts = deserialize(joinpath(OUT, "starts.jls"))

# the quantities of a state of a case (W with the factor-7 envelope at the T of the state)
function record(case, u)
    M = moments(case, u)
    (T = vals(case, u, :cold_electron_temperature), n_cold = vals(case, u, :cold_electron_density),
        n_below = M.Nbelow, n_hot = vals(case, u, :hot_electron_density), j_hot = vals(case, u, :hot_electron_current_density),
        j_ohm = vals(case, u, :ohmic_current_density), E = vals(case, u, :electric_field),
        W_below = M.Wbelow, W_hot = M.Whot, W_cold = vals(case, u, :cold_thermal_energy), W_all = M.Wall, n_f = M.Nall)
end

results = Dict{Symbol, Any}()
function dorun(key, case, u0, label)
    println("\n================ run $key: $label ================"); flush(stdout)
    t0 = time()
    system, sol, msg = run_case(case, u0, times; label = string(key))
    if sol === nothing
        results[key] = (; label, failed = msg)
        return nothing
    end
    println("  wall time ", round(time() - t0; digits = 1), " s; ", stats_line(sol))
    recs = [record(case, sol.states[:, k]) for k in 1:length(sol.times)]
    results[key] = (; label, case, sol, recs, times = sol.times)
    return sol
end

# budget of the energy of f per step for a prototype run
function budget(pcase, sol)
    println("  energy budget of f per step [W/m^3] (whole-grid energy moment; backward Euler: terms at u_new, cooling at u_old)")
    layout = pcase.layout
    fp = fp_of(pcase)
    feq = eq_of(pcase, :f)
    rows = []
    for k in 2:length(sol.times)
        u_old, u, t = sol.states[:, k - 1], sol.states[:, k], sol.times[k]
        dt = t - sol.times[k - 1]
        prepare_terms!(pcase, u_old, t)
        dW = (moments(pcase, u).Wall .- moments(pcase, u_old).Wall) ./ dt
        dN = (moments(pcase, u).Nall .- moments(pcase, u_old).Nall) ./ dt
        parts = Dict{String, Vector{Float64}}(); dens = Dict{String, Vector{Float64}}()
        local c
        for term in feq.terms
            name = string(nameof(typeof(term)))
            if term isa Cooling
                c = term; c.key = nothing; compute!(c, u_old, layout, t)
                W, N = dmoments(pcase, (o = zeros(length(u)); for ir in 1:pcase.nr
                    kinetic_view(o, layout, unknown_id(layout, :f), ir) .+= c.df[:, :, ir]
                end; o))
            else
                W, N = term_f(pcase, term, u, t)
            end
            parts[name] = W; dens[name] = N
        end
        Lfe = tv(pcase, KineticEnergyLossTerm(fp; part = :free_electrons, coefficient = 1.0), u, t)
        Ls = tv(pcase, KineticEnergyLossTerm(fp; part = :sinks, coefficient = 1.0), u, t)
        out = tv(pcase, MomentumBoundaryFluxTerm(:f, fp; weight = :energy), u, t)
        field = parts["FokkerPlanckTerm"] .+ Lfe .+ Ls .+ out
        T = vals(pcase, u, :cold_electron_temperature); n = vals(pcase, u, :cold_electron_density)
        uM = withf(pcase, u, maxwellian_array(pcase, n, T))
        LfeM = tv(pcase, KineticEnergyLossTerm(fp; part = :free_electrons, coefficient = 1.0), uM, t)
        sumterms = sum(values(parts))
        println(@sprintf("   step %d t=%.2e dt=%.1e", k - 1, t, dt))
        pr("    dW_f/dt", dW)
        pr("    sum of the term moments", sumterms)
        pr("      field work (FP − collision parts)", field)
        pr("      e-e collisions on the free electrons (−L_fe)", -Lfe)
        pr("      sinks (−L_s)", -Ls)
        pr("      outflow at p_max (−)", -out)
        pr("      cooling (P_fluid)", parts["Cooling"])
        pr("      particle source", parts["ParticleSourceTerm"])
        pr("      knock-on source", parts["RosenbluthPutvinskiSourceTerm"])
        for kk in sort(collect(keys(c.parts)))
            pr("        P_fluid: $kk", c.parts[kk])
        end
        pr("    −L_fe on a Maxwellian at (n_cold,T) of u_new", -LfeM)
        pr("    cooling a [1/s], CFL a p g J/(JΔp) dt", c.a); pr("    CFL", c.rate .* dt)
        pr("    density moment of the cooling [m^-3/s]", dens["Cooling"])
        pr("    dn_f/dt vs particle source density", [dN; dens["ParticleSourceTerm"]])
        push!(rows, (; dW, sumterms, field, ee = -Lfe, sinks = -Ls, out = -out, cool = parts["Cooling"],
            src = parts["ParticleSourceTerm"], ko = parts["RosenbluthPutvinskiSourceTerm"], eeM = -LfeM, P = copy(c.P)))
    end
    flush(stdout)
    return rows
end

# bulk shape at the end: f_iso / Maxwellian(n_cold, T) below the threshold
function shape_report(case, u, label)
    T = vals(case, u, :cold_electron_temperature); n = vals(case, u, :cold_electron_density)
    F = farr(case, u); FM = maxwellian_array(case, n, T)
    region = region_of(case)
    println("  bulk shape at the end, $label: f_iso/f_M per p cell with p < p0, and the L1 distance Σ(1−w)|f_iso−f_M| / n")
    for ir in 1:case.nr
        m = momentum_grid(case.grid, ir); dxi = m.xi_widths; D, _ = weights(case.grid, ir)
        w = envelope(region, case.grid, ir, T[ir])
        iso(A) = [sum(A[i, :, ir] .* dxi) / sum(dxi) for i in 1:momentum_cell_count(m)]
        fi, mi = iso(F), iso(FM)
        Dp = vec(sum(D; dims = 2))
        L1 = sum((1 .- w) .* Dp .* abs.(fi .- mi)) / n[ir]
        p0 = p0_of(region, T[ir])
        sel = [i for i in eachindex(fi) if m.p_centers[i] < p0]
        println(@sprintf("   r%d T=%.4g eV p0=%.3g L1=%.3e : ", ir, T[ir], p0, L1),
            join([@sprintf("%.3g:%.3g", m.p_centers[i], fi[i] / mi[i]) for i in sel], " "))
    end
    flush(stdout)
end

budgets_ = Dict{Symbol, Any}()
if :pa in VARIANTS
    sol = dorun(:pa, pc, starts.ua, "(i) prototype, start (a): T from the moment of DREAM's f")
    sol === nothing || (budgets_[:pa] = budget(pc, sol); shape_report(pc, sol.states[:, end], "(i)"))
end
if :pb in VARIANTS
    sol = dorun(:pb, pc, starts.ub, "(ii) prototype, start (b): Maxwellian bulk at DREAM's (n_cold, T_cold)")
    sol === nothing || (budgets_[:pb] = budget(pc, sol); shape_report(pc, sol.states[:, end], "(ii)"))
end
if :pas in VARIANTS
    sol = dorun(:pas, pcs, starts.ub, "(ii') prototype, start (b), cooling compensates the sinks below the threshold")
    sol === nothing || (budgets_[:pas] = budget(pcs, sol); shape_report(pcs, sol.states[:, end], "(ii')"))
end
if :cold in VARIANTS
    sol = dorun(:cold, base_r, initial_state(base_r, START), "(iii) runes2p, cold-energy row")
    sol === nothing || shape_report(base_r, sol.states[:, end], "(iii)")
end
if :total in VARIANTS
    tc = total_case(base_r)
    dorun(:total, tc, initial_state(tc, START), "(iv) runes2p, total-energy row")
end
if :dream in VARIANTS
    bd = base_case(:dream)
    dorun(:dream, bd, initial_state(bd, START), "(v) dream, cold-energy row")
end

# DREAM's values (vi), W with the factor-7 envelope at DREAM's T
function dream_records()
    nt = min(length(TIMES), length(times))
    recs = []
    for k in 1:nt
        T = START["T_cold"][:, k]
        F = START["f_hot"][:, :, :, k]
        M = fmoments(base_r, F, T)
        push!(recs, (T, n_cold = START["n_cold"][:, k], n_below = M.Nbelow, n_hot = START["n_hot"][:, k], j_hot = START["j_hot"][:, k],
            j_ohm = START["j_ohm"][:, k], E = START["E_field"][:, k], W_below = M.Wbelow, W_hot = M.Whot,
            W_cold = START["W_cold"][:, k], W_all = M.Wall, n_f = M.Nall))
    end
    recs
end
results[:file] = (; label = "(vi) DREAM's file", recs = dream_records(), times = TIMES[1:length(dream_records())])

# tables
order = [k for k in (:pa, :pb, :pas, :cold, :total, :dream, :file) if haskey(results, k) && !haskey(results[k], :failed)]
println("\n================ tables (per radius r1 | r2) ================")
println("runs: ", join(["$k = $(results[k].label)" for k in order], "; "))
units = Dict(:T => "eV", :n_cold => "m^-3", :n_below => "m^-3", :n_hot => "m^-3", :j_hot => "A/m^2", :j_ohm => "A/m^2", :E => "V/m",
    :W_below => "J/m^3", :W_hot => "J/m^3", :W_cold => "J/m^3", :W_all => "J/m^3", :n_f => "m^-3")
for q in (:T, :n_cold, :n_below, :n_hot, :j_hot, :j_ohm, :E, :W_below, :W_hot, :W_cold)
    println("\n", q, " [", units[q], "]   columns: ", join(["$k" for k in order], " | "))
    for k in 1:length(times)
        cells = String[]
        for key in order
            r = results[key]
            push!(cells, k <= length(r.recs) ? join([@sprintf("%.4e", x) for x in getfield(r.recs[k], q)], " ") : "      -")
        end
        println(@sprintf("  t=%.1e  ", times[k]), join(cells, " | "))
    end
end
for (ref, name) in ((:file, "DREAM (vi)"), (:cold, "(iii) cold row"))
    haskey(results, ref) || continue
    println("\nrelative difference of T_cold to $name, (T − T_ref)/T_ref   columns: ", join(["$k" for k in order], " | "))
    for k in 1:length(times)
        k <= length(results[ref].recs) || continue
        Tr = results[ref].recs[k].T
        cells = [k <= length(results[key].recs) ? join([@sprintf("%+.3e", x) for x in (results[key].recs[k].T .- Tr) ./ Tr], " ") : "      -" for key in order]
        println(@sprintf("  t=%.1e  ", times[k]), join(cells, " | "))
    end
end
for k in order
    haskey(results[k], :sol) || continue
    s = results[k].sol.statistics
    println("\nNewton/steps $k: ", stats_line(results[k].sol))
end
for k in keys(results)
    haskey(results[k], :failed) && println("\nFAILED $k: ", results[k].failed)
end
serialize(joinpath(OUT, "results_$(get(ENV, "PROTO_SCHEME", "tcdf"))_$(NSTEPS)_$(join(VARIANTS, "-")).jls"),
    Dict(k => (haskey(v, :recs) ? (v.label, v.recs, v.times) : (v.label, v.failed)) for (k, v) in results))
serialize(joinpath(OUT, "budgets_$(get(ENV, "PROTO_SCHEME", "tcdf"))_$(NSTEPS)_$(join(VARIANTS, "-")).jls"), budgets_)
println("END")
