# The disruption case with DREAM's runaway grid above the hot-tail grid (`build` of
# example_disruption.jl with `runaway = true`; the case script's --runaway-grid,
# reference/dream/cases/disruption_kinetic.py, with its default grid, the constants `RUNAWAY_*`
# of example_disruption.jl). The system starts from DREAM's state at the start of the quench of
# dream_disruption_small.h5, a run without the runaway grid, with f_re zero, and takes 3 of its
# steps of 1e-7 s at both compat values. f_re grows from zero to 1.8e7 in them, so its scale is
# `STARTUP_SCALE` rather than the default near DREAM's f_re in the quench (1e16): Newton's tests
# on f_re's rows are then at 1e-10 of it, and a step that left f_re unchanged would fail them.
# Nothing is compared with DREAM here. The checks, independent of DREAM: the steps are accepted,
# the budgets of the case hold within the bounds of the phases started from DREAM's state
# (`budget_checks`, among them ∫f − n_cold − n_hot fixed to round-off), the electrons that the
# interface takes from f are those it gives to f_re and to the row of n_re, and n_re − ∫f_re grows
# by the avalanche beyond f_re's p_max and the outflow there (`runaway_budgets`). Without a runaway
# grid the builder gives the unknowns and the channels of the row of W_cold of the case as before.
# The start state takes the particle source S_particle of a DREAM output that holds it
# (`initial_state`). Needs the ADF11 cache.
using HDF5
include(joinpath(@__DIR__, "example_disruption.jl"))

# the scale of f_re in these steps, near its largest value after them
const STARTUP_SCALE = 1e7

# the bounds of `runaway_budgets`, at least 45 times the values measured at -O0 (in the comment;
# compat = dream / runes2p).
# The particles through the interface agree to round-off. The balance of n_re − ∫f_re is limited
# by the change of n_re in a step, 3e-6 of n_re at the first radius: the round-off of n_re,
# 1e-16 of it, is then 4e-11 of the change
#! format: off
const RUNAWAY_BUDGETS = (interface = 2e-14,   # 2.0e-16 / 1.8e-16
                         avalanche = 1e-9)    # 2.2e-11 / 2.1e-11
#! format: on

@testset "the disruption case with a runaway grid" begin
    if !isdir(adas_cache_directory())
        @info "no ADF11 cache at $(adas_cache_directory()): the disruption case with a runaway grid is skipped"
        @test_skip isdir(adas_cache_directory())
    else
        file = reference_file("dream_disruption_small.h5")
        start = h5open(h5 -> read_dream_state(h5["quench"]), file, "r")
        times = collect(start["t"])[1:4]
        unknown_names(case) = [unknown.name for unknown in case.layout.unknowns]

        @testset "without a runaway grid" begin
            # a file without DREAM's runaway grid builds the case without one, as before
            @test !haskey(start, "p1_f_re")
            plain = build(start)
            explicit = build(start; runaway = false)
            @test plain.runaway_grid === nothing && plain.interface === nothing
            @test unknown_names(plain) == unknown_names(explicit)
            @test :f_re ∉ unknown_names(plain)
            @test plain.channel_labels == ["E j_ohm", "radiation and binding",
                "heat diffusion", "e c E_c n_re", "exchange (electrons)"]
            # the runaway grid adds f_re alone to the unknowns
            with_grid = build(start; runaway = true)
            @test unknown_names(with_grid) == vcat(unknown_names(plain), :f_re)
        end

        @testset "the particle source of the start" begin
            # the comparison files of the reference sets hold no S_particle, so the start state
            # takes a source of distinct values per radius and time, written as DREAM's output
            # holds it and read back by read_dream_state
            @test !haskey(start, "S_particle")
            case = build(start; runaway = true)
            source_id = unknown_id(case.layout, :particle_source)
            without = initial_state(case, start)
            @test all(iszero, unknown_values(without, case.layout, source_id))
            source = [1e24 * (ir + 0.1 * it)
                      for ir in 1:(case.nr), it in eachindex(start["t"])]
            file = tempname() * ".h5"
            h5open(reference_file("dream_disruption_small.h5"), "r") do h5
                h5open(file, "w") do out
                    copy_object(h5["quench"], out, "quench")
                    out["quench/eqsys/S_particle"] = source
                end
            end
            with_source = h5open(h5 -> read_dream_state(h5["quench"]), file, "r")
            rm(file)
            @test with_source["S_particle"] == source
            with = initial_state(case, with_source)
            @test unknown_values(with, case.layout, source_id) == source[:, 1]
            others = trues(length(with))
            for ir in 1:(case.nr)
                others[unknown_range(case.layout, source_id, ir)] .= false
            end
            @test with[others] == without[others]
            # the superthermal mode has no particle source and ignores S_particle
            superthermal = build(with_source; mode = :superthermal, runaway = true)
            @test :particle_source ∉ unknown_names(superthermal)
            @test initial_state(superthermal, with_source) ==
                  initial_state(superthermal, start)
        end

        @testset "compat = $compat" for compat in (:dream, :runes2p)
            case = build(start; compat, runaway = true, runaway_scale = STARTUP_SCALE)
            p_faces = momentum_grid(case.runaway_grid, 1).p_faces
            @test first(p_faces) == last(start["p1_f"]) && last(p_faces) == RUNAWAY_P_MAX
            @test length(p_faces) == 33
            @test pitch_cell_count(momentum_grid(case.runaway_grid, 1)) ==
                  RUNAWAY_PITCH_CELLS
            @test case.interface.pitch_mapping === (compat === :dream ? :dream : :overlap)
            # the scheme of f_re: that of f at compat = dream, exponential fitting at runes2p,
            # where the flux through the interface is upwind next to the four-point scheme of f
            @test case.fp_re.operator.p_scheme ===
                  (compat === :dream ? :tcdf : :exponential_fitting)
            if compat === :runes2p
                @test case.interface.scheme === :upwind
            end
            u0 = initial_state(case, start)
            re_id = unknown_id(case.layout, :f_re)
            @test all(iszero, kinetic_array(u0, case.layout, re_id))
            system = EquationSystem(case.layout, case.equations;
                time_scale = times[2] - times[1], initial_state = u0)
            solution = integrate(system, u0, times; stepper = FixedSteps())
            @test all(solution.statistics.accepted)
            label = "quench with a runaway grid, compat = $compat"
            for (_, value, bound) in budget_checks(case, solution, label, DRIFT["quench"])
                @test value < bound
            end
            # the closures of the runaway grid, and electrons cross into f_re
            for (_, value, bound) in runaway_budget_checks(case, solution, label, RUNAWAY_BUDGETS)
                @test value < bound
            end
            # f_re fills
            final = density(case.runaway_grid, kinetic_array(solution.states[:, end], case.layout, re_id))
            @test all(>(0), final)
        end
    end
end
