# Prototype for issue #312: T_cold from the energy moment of f below the threshold of the hot region
# (fully kinetic mode, compat = :runes2p), with the fluid channels of the cold row applied to f
# as a shape-preserving cooling advection. Scratch code, outside the package.
# the checkout that holds this folder, and the folder of the serialized results
const WT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUT = mkpath(joinpath(@__DIR__, "output"))
const TESTDIR = joinpath(WT, "test")
using Runes2P, HDF5, Printf, LinearAlgebra
include(joinpath(TESTDIR, "helpers.jl")); using .TestHelpers
include(joinpath(TESTDIR, "example_cases.jl"))
include(joinpath(TESTDIR, "example_disruption.jl"))
println("HEAD ", readchomp(`git -C $WT rev-parse HEAD`))
println("pathof ", pathof(Runes2P)); flush(stdout)

const e = Runes2P.Constants.ELEMENTARY_CHARGE
const MC2 = Runes2P.Constants.ELECTRON_REST_ENERGY   # eV
const FILE = reference_file("dream_disruption_small.h5")
const START = h5open(h5 -> read_dream_state(h5["quench"]), FILE, "r")
const TIMES = collect(START["t"])

# ---------------------------------------------------------------- helpers on a case
eq_of(case, name) = only(e for e in case.equations if e.unknown === name)
fp_of(case) = only(t for t in eq_of(case, :f).terms if t isa FokkerPlanckTerm)
rp_of(case) = only(t for t in eq_of(case, :f).terms if t isa RosenbluthPutvinskiSourceTerm)
region_of(case) = only(t for t in eq_of(case, :hot_electron_density).terms if t isa DensityMomentTerm).envelopes.region
vals(case, u, name) = collect(unknown_values(u, case.layout, unknown_id(case.layout, name)))
setvals!(case, u, name, v) = (unknown_values(u, case.layout, unknown_id(case.layout, name)) .= v; u)
farr(case, u) = copy(kinetic_array(u, case.layout, unknown_id(case.layout, :f)))
withf(case, u, f) = (s = copy(u); kinetic_array(s, case.layout, unknown_id(case.layout, :f)) .= f; s)
v(x) = join([@sprintf("%.4e", y) for y in x], ", ")
pr(label, x) = (println(@sprintf("  %-60s [%s]", label, v(x))); flush(stdout))
p0_of(region, T) = region.limit.factor * sqrt(2T / MC2)

# energy weight ε J Δp Δξ / V' and density weight per cell at radius ir
function weights(grid, ir)
    m = momentum_grid(grid, ir); J = phase_space_jacobian(grid, ir); Vr = grid.radial_grid.spatial_jacobian[ir]
    np, nx = momentum_cell_count(m), pitch_cell_count(m)
    D = [J[i, j] * m.p_widths[i] * m.xi_widths[j] / Vr for i in 1:np, j in 1:nx]
    E = [Runes2P._kinetic_energy(m.p_centers[i]) for i in 1:np] .* D
    return D, E
end
envelope(region, grid, ir, T) = Runes2P.momentum_envelope!(zeros(momentum_cell_count(momentum_grid(grid, ir))), region,
    momentum_grid(grid, ir), p0_of(region, T))
# moments of an f array (np, nξ, nr) with the hot envelope at T per radius
function fmoments(case, F, T)
    g = case.grid; region = region_of(case); nr = case.nr
    out = (Wall = zeros(nr), Whot = zeros(nr), Nall = zeros(nr), Nhot = zeros(nr))
    for ir in 1:nr
        D, E = weights(g, ir); w = envelope(region, g, ir, T[ir])
        f = F[:, :, ir]
        out.Wall[ir] = sum(E .* f); out.Whot[ir] = sum(w .* E .* f)
        out.Nall[ir] = sum(D .* f); out.Nhot[ir] = sum(w .* D .* f)
    end
    return merge(out, (Wbelow = out.Wall .- out.Whot, Nbelow = out.Nall .- out.Nhot))
end
moments(case, u) = fmoments(case, farr(case, u), vals(case, u, :cold_electron_temperature))
# energy and density moment (whole grid) of a change df on the f row
function dmoments(case, out)
    F = kinetic_array(out, case.layout, unknown_id(case.layout, :f)); g = case.grid
    W = zeros(case.nr); N = zeros(case.nr)
    for ir in 1:case.nr
        D, E = weights(g, ir); W[ir] = sum(E .* F[:, :, ir]); N[ir] = sum(D .* F[:, :, ir])
    end
    return W, N
end
term_f(case, term, u, t) = dmoments(case, Runes2P.rhs!(zeros(length(u)), term, u, case.layout, unknown_id(case.layout, :f), t))
radial_row(case) = unknown_id(case.layout, :cold_thermal_energy)
tv(case, term, u, t) = _term_values(term, u, case.layout, radial_row(case), t)

# ---------------------------------------------------------------- the cooling term
# Explicit term on f: the advection A_p = a_k p g(p) in momentum, g = 1 − w at the faces (zero in
# the hot region), donor cell by the sign of a_k, a_k per radius chosen so that the energy moment
# of the term is P_fluid,k evaluated at the state of the start of the step (the state the system
# passes to an explicit term). Flux form: no flux at p = 0 and where g = 0, so it moves no electrons.
mutable struct Cooling <: Runes2P.EquationTerm
    grid::Any
    region::Any                  # the hot region (above), for the envelope at T_old
    channels::Vector{Pair{String, Any}}   # the fluid channels of the cold row (W/m^3)
    fp::Any                      # the Fokker-Planck term of f, for the optional sinks compensation
    compensate_sinks::Bool
    key::Any
    df::Array{Float64, 3}
    a::Vector{Float64}
    P::Vector{Float64}
    parts::Dict{String, Vector{Float64}}
    rate::Vector{Float64}        # max over cells of |a| p_f g J_f / (J_c Δp) [1/s]
end
function Cooling(case; compensate_sinks = false)
    T_eq = eq_of(case, :cold_electron_temperature)
    labels = case.channel_labels
    ch = Pair{String, Any}[l => t for (l, t) in zip(labels, T_eq.terms) if l != "E j_ohm"]
    g = case.grid; np, nx = momentum_cell_count(momentum_grid(g, 1)), pitch_cell_count(momentum_grid(g, 1))
    Cooling(g, region_of(case), ch, fp_of(case), compensate_sinks, nothing, zeros(np, nx, case.nr),
        zeros(case.nr), zeros(case.nr), Dict{String, Vector{Float64}}(), zeros(case.nr))
end
Runes2P.dependencies(::Cooling) = [:f]
Runes2P.treatment(::Cooling) = :explicit
Runes2P.sparsity!(entries, ::Cooling, layout, row) = entries

# the unit advection (|a| = 1) with sign s on f at radius ir, as df (np, nξ)
function unit_advection(grid, ir, f, w, s)
    m = momentum_grid(grid, ir); Jc = phase_space_jacobian(grid, ir); Jf = phase_space_jacobian_p_faces(grid, ir)
    np, nx = size(f)
    pf = m.p_faces
    gf = [k == 1 ? 1.0 : k == np + 1 ? 0.0 : 1 - (w[k - 1] + w[k]) / 2 for k in 1:(np + 1)]
    Φ = zeros(np + 1, nx)     # flux through face k (lower face of cell k), times J_f
    for j in 1:nx, k in 2:np
        donor = s > 0 ? f[k - 1, j] : f[k, j]
        Φ[k, j] = s * pf[k] * gf[k] * Jf[k, j] * donor
    end
    df = [-(Φ[i + 1, j] - Φ[i, j]) / (Jc[i, j] * m.p_widths[i]) for i in 1:np, j in 1:nx]
    rate = maximum(pf[k] * gf[k] * Jf[k, j] / (Jc[kk, j] * m.p_widths[kk]) for j in 1:nx, k in 2:np for kk in (k - 1, k))
    return df, rate
end

function compute!(c::Cooling, u, layout, t)
    key = (t, hash(u))
    key == c.key && return c
    Trow = unknown_id(layout, :cold_thermal_energy)
    nr = radial_cell_count(layout)
    parts = Dict{String, Vector{Float64}}()
    for (label, term) in c.channels
        Runes2P.prepare_step!(term, u, layout, unknown_id(layout, :cold_electron_temperature), t)
        parts[label] = _term_values(term, u, layout, Trow, t)
    end
    E = collect(unknown_values(u, layout, unknown_id(layout, :electric_field)))
    johm = collect(unknown_values(u, layout, unknown_id(layout, :ohmic_current_density)))
    jhot = collect(unknown_values(u, layout, unknown_id(layout, :hot_electron_current_density)))
    jf = _term_values(CurrentMomentTerm(:f, c.grid; pitch_rule = :cell_center), u, layout, Trow, t)
    parts["E^2(sigma-sigma_kin)"] = E .* (johm .- (jf .- jhot))
    if c.compensate_sinks
        below = MomentumRegion(c.region.limit; side = :below, envelope = c.region.envelope)
        parts["+sinks below (compensation)"] = _term_values(
            KineticEnergyLossTerm(c.fp; part = :sinks, region = below, coefficient = 1.0), u, layout, Trow, t)
    end
    P = sum(values(parts))
    T = collect(unknown_values(u, layout, unknown_id(layout, :cold_electron_temperature)))
    fid = unknown_id(layout, :f)
    for ir in 1:nr
        f = Array(kinetic_view(u, layout, fid, ir))
        w = envelope(c.region, c.grid, ir, T[ir])
        s = P[ir] >= 0 ? 1.0 : -1.0
        df, rate = unit_advection(c.grid, ir, f, w, s)
        _, Ew = weights(c.grid, ir)
        Eunit = sum(Ew .* df)
        a = P[ir] / Eunit        # >= 0 when the unit advection has the sign of P
        c.df[:, :, ir] .= a .* df
        c.a[ir] = s * a; c.rate[ir] = abs(a) * rate
    end
    c.P .= P; c.parts = parts; c.key = key
    return c
end
function Runes2P.rhs!(out, c::Cooling, u, layout, row, t)
    compute!(c, u, layout, t)
    for ir in 1:radial_cell_count(layout)
        kinetic_view(out, layout, row, ir) .+= c.df[:, :, ir]
    end
    return out
end

# ---------------------------------------------------------------- case variants
base_case(compat) = build(START; compat, runaway = false, scheme = Symbol(get(ENV, "PROTO_SCHEME", "tcdf")))

function below_region(case)
    r = region_of(case)
    MomentumRegion(r.limit; side = :below, envelope = r.envelope, smooth_width = r.smooth_width)
end

# the prototype: T from the energy moment below the threshold, the cooling on f
function proto_case(case; compensate_sinks = false)
    below = below_region(case)
    cooling = Cooling(case; compensate_sinks)
    Trow = Equation(:cold_electron_temperature, Algebraic(),
        EnergyMomentTerm(:f, case.grid; region = below),
        ProductTerm(:cold_electron_density, :cold_electron_temperature; coefficient = -THERMAL_ENERGY_FACTOR);
        scale = 1e4)
    feq = eq_of(case, :f)
    equations = Equation[]
    for eq in case.equations
        if eq.unknown === :cold_electron_temperature
            push!(equations, Trow)
        elseif eq.unknown === :f
            push!(equations, Equation(:f, TimeDerivative(), feq.terms..., cooling))
        else
            push!(equations, eq)
        end
    end
    # the electron side of the ion exchange left the row of T; the system requires both sides, so
    # it goes into a diagnostic unknown x = Q_e (W/m^3) that nothing reads; the cooling applies Q_e to f
    exch = only(t for (l, t) in zip(case.channel_labels, eq_of(case, :cold_electron_temperature).terms) if l == "exchange (electrons)")
    push!(equations, Equation(:exchange_electrons, Algebraic(), exch, LinearTerm(:exchange_electrons; coefficient = -1)))
    layout = UnknownLayout(case.radial, vcat(case.layout.unknowns, RadialUnknown(:exchange_electrons; scale = 1e6)))
    return merge(case, (; layout, equations, cooling, below))
end

# the case on the row of the total energy (total_energy_equation), for the comparison
function total_case(case)
    region = region_of(case); fp, rp = fp_of(case), rp_of(case)
    cold_row = eq_of(case, :cold_electron_temperature)
    row_terms = [collect(cold_row.terms); ProductTerm(:electric_field, :hot_electron_current_density);
                 KineticEnergyLossTerm(fp; part = :sinks, region);
                 MomentumBoundaryFluxTerm(:f, fp; weight = :energy, coefficient = -1);
                 KineticEnergyChangeTerm(rp)]
    unknowns = vcat(case.layout.unknowns, RadialUnknown(:total_electron_energy; scale = 1e4))
    layout = UnknownLayout(case.radial, unknowns)
    equations = [e for e in case.equations if e.unknown !== :cold_electron_temperature]
    push!(equations, total_energy_equation(:f, case.grid; region))
    push!(equations, Equation(:cold_electron_temperature, TimeDerivative(:total_electron_energy), row_terms...))
    return merge(case, (; layout, equations))
end

# solve the linear algebraic rows of `names` for their unknowns, entry by entry, at the state u
function solve_rows!(case, u, names; t = TIMES[1], passes = 3)
    layout = case.layout
    for _ in 1:passes, name in names
        eq = eq_of(case, name)
        id = unknown_id(layout, name)
        rowvals_(uu) = (out = zeros(length(uu)); for term in eq.terms
            Runes2P.rhs!(out, term, uu, layout, id, t)
        end; collect(unknown_values(out, layout, id)))
        x = unknown_values(u, layout, id)
        x0 = copy(x); F0 = rowvals_(u)
        s = layout.unknowns[id].scale
        x .= x0 .+ s; F1 = rowvals_(u)
        c = (F1 .- F0) ./ s
        x .= x0 .- F0 ./ c
    end
    return u
end
const DERIVED = [:hot_electron_density, :cold_electron_density, :hot_electron_current_density,
    :ohmic_current_density, :total_current_density, :cold_thermal_energy]

# the residual of the moment row per radius at T (rows solved), relative to 3/2 e n T
function Frow(pc, u, T)
    s = setvals!(pc, copy(u), :cold_electron_temperature, T)
    solve_rows!(pc, s, DERIVED)
    M = moments(pc, s); n = vals(pc, s, :cold_electron_density)
    return M.Wbelow .- THERMAL_ENERGY_FACTOR .* n .* T, s
end
# per-radius bisection of the moment row in T between lo and hi
function moment_temperature(pc, u, lo, hi; iters = 80)
    flo, _ = Frow(pc, u, lo); fhi, _ = Frow(pc, u, hi)
    all(sign.(flo) .!= sign.(fhi)) || error("no bracket: F(lo)=$flo F(hi)=$fhi")
    lo, hi = copy(lo), copy(hi)
    for _ in 1:iters
        mid = (lo .+ hi) ./ 2
        fm, _ = Frow(pc, u, mid)
        for ir in eachindex(mid)
            if sign(fm[ir]) == sign(flo[ir])
                lo[ir] = mid[ir]; flo[ir] = fm[ir]
            else
                hi[ir] = mid[ir]
            end
        end
    end
    T = (lo .+ hi) ./ 2
    _, s = Frow(pc, u, T)
    return T, s
end

# Maxwellian array at (n, T) per radius on the grid of the case
function maxwellian_array(case, n, T)
    u = zeros(length(case.layout)); fill_maxwellian!(u, case.layout, :f, n, T); farr(case, u)
end

function run_case(case, u0, times; stepper = FixedSteps(), label = "")
    system = EquationSystem(case.layout, case.equations; time_scale = times[2] - times[1], initial_state = u0)
    u = copy(u0)
    try
        sol = integrate(system, u, times; stepper)
        return system, sol, nothing
    catch err
        msg = first(sprint(showerror, err), 1500)
        println("RUN FAILED [$label]: ", msg); flush(stdout)
        return system, nothing, msg
    end
end
function stats_line(sol)
    s = sol.statistics
    "steps attempted=$(length(s.accepted)) accepted=$(count(s.accepted)) not accepted=$(count(!, s.accepted)) " *
    "newton=$(s.newton_iterations) status=$(unique(s.status))"
end
