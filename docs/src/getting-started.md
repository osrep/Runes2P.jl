# Getting started

```@meta
CurrentModule = Runes2P
```

This page covers the installation, a first run with the examples, the settings of a run and its output files.

## Installation

Runes2P needs Julia 1.13.
It is not registered in the General registry and is added from its repository:

```julia
using Pkg
Pkg.add(url = "https://github.com/osrep/Runes2P.jl")
```

The HDF5 output is a package extension.
It loads with `using HDF5` next to `using Runes2P` and defines [`write_output`](@ref) and [`read_output`](@ref).
Checkpoints need it loaded.

The examples run from a checkout of the repository.
Their environment, `examples/Project.toml`, takes Runes2P from the parent directory and adds HDF5:

```sh
julia --project=examples -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=examples examples/relaxation.jl
julia --project=examples examples/cooling.jl scheme=van_albada solver=gmres
julia --project=examples examples/avalanche.jl treatment=implicit solver=gmres
```

The scripts are described in [The examples](#The-examples), and the arguments such as `scheme=` in [Switches](#Switches).

## The examples

Each script solves one kinetic equation for the electron distribution ``f(p, \xi, r)`` and runs without arguments.
Its header lists its default values and further command lines.

| Script | What it integrates | What it reports |
|---|---|---|
| `relaxation.jl` | a Maxwell-Jüttner distribution in a fixed plasma and field | with `p_star=`, the fraction of the density above that momentum |
| `cooling.jl` | a thermal quench that starts in the slide-away regime | the saves at which ``f`` has negative cells, and the first of them in momentum |
| `avalanche.jl` | knock-on multiplication of a Maxwell-Jüttner seed | the mean growth rate of the density over the run |

### `relaxation.jl`

- The initial state is a Maxwell-Jüttner distribution at the plasma temperature, or at `T0=`.
- The plasma is deuterium, and the plasma and the field stay fixed.
- The operator is that of [`fokker_planck_operator`](@ref).
- `p_star=` adds the fraction of the density above that momentum at every save.
- The plasma values are choices of this package ([Choices of this package](choices.md)).

### `cooling.jl`

- The temperature decays exponentially from `T=` to `T_final=` over `t_end=`.
- [`prescribed_fokker_planck_term`](@ref) rebuilds the collision coefficients from the temperature.
- The field stays at `e_over_ec=` times the critical field at the initial temperature.
- With `nr=` above one, radial diffusion along stochastic field lines couples the radii ([`add_rechester_rosenbluth_diffusion!`](@ref)).
  Its strength is `dbb=` (``\delta B/B``), and ``f = 0`` at the wall.
  The diffusion is that of A. B. Rechester and M. N. Rosenbluth, Phys. Rev. Lett. 40, 38 (1978), doi:10.1103/PhysRevLett.40.38.
- The plasma is that of `reference/dream/cases/cooling_positivity.py`, cooled over 2 ms where that case takes 1 ms.
- The radial values are those of DREAM's SPI example.
- The grid and the time step are a choice of this package.

### `avalanche.jl`

- The seed is a Maxwell-Jüttner distribution of temperature 100 keV.
- The knock-on collisions are those of [`KnockOnOperator`](@ref) and [`KnockOnTerm`](@ref).
- The setup is that of `examples/avalanche/generateCylindrical.py` of DREAM's `boltzmann-eqterm` branch.
- The script takes the values of that example, among them these.
  - Superthermal collisions.
  - The energy-dependent Coulomb logarithm.
  - DREAM's ion-size table.
- The interpolation scheme is the default of this package, where the example uses TCDF (`scheme=tcdf`).

The switch `treatment=` of [`KnockOnTerm`](@ref) sets the time level of the knock-on term.

| `treatment=` | The knock-on term |
|---|---|
| `explicit_primary` (default) | the primary distribution at the previous time level, as in DREAM |
| `explicit` | the whole term at the previous time level |
| `implicit` | the whole term at the new time level |

`implicit` applies the term matrix-free in each Krylov iteration of [`GMRES`](@ref).
The direct solvers do not take it, so with `treatment=implicit` the solver defaults to GMRES.

### What every script prints

Every script prints the following.

- The segment and its status.
- The solver statistics ([`SolverStatistics`](@ref)).
- The negative entries of the saved states.
- The density per radius at every save.

A segment is one part of a run that is split over several allocations ([Checkpoints, segments and the stop](model/time-stepping.md#Checkpoints,-segments-and-the-stop)).
In a GMRES run the statistics include the Krylov iterations and factorizations.
Before the timed run a script integrates over the first three time points, so that the timing leaves out the compilation.
`relaxation.jl` and `cooling.jl` end with [`exit_status`](@ref) as the exit status of the process.

## The parts of a run

The three scripts build a run in the same order.

1. A [`MomentumGrid`](@ref) from momentum faces ([`uniform_p_faces`](@ref), [`biuniform_p_faces`](@ref), [`tanh_p_faces`](@ref)) and pitch faces ([`uniform_xi_faces`](@ref), [`biuniform_theta_xi_faces`](@ref)), combined with a [`RadialGrid`](@ref) by [`cylindrical_grid`](@ref) ([Grids](model/grids.md)).
2. A [`PlasmaState`](@ref) with its [`IonSpecies`](@ref), and the [`CollisionSettings`](@ref) of the collision frequencies ([Collision models](model/collisions.md)).
3. The kinetic operator, an [`AdvectionDiffusion`](@ref) with its interpolation schemes and boundaries ([Positivity and flux limiters](model/positivity.md)).
4. An [`UnknownLayout`](@ref) with one [`KineticUnknown`](@ref) `:f`, and an [`EquationSystem`](@ref) of one [`Equation`](@ref): a [`TimeDerivative`](@ref) and the terms, among them [`FokkerPlanckTerm`](@ref) or [`prescribed_fokker_planck_term`](@ref) ([Unknowns and the equation system](model/unknowns.md)).
5. A [`NewtonSolver`](@ref) with its linear solver, and [`integrate`](@ref) over the time points, which returns a [`TimeSolution`](@ref) ([Time stepping and solvers](model/time-stepping.md)).

The shared helpers of the scripts are in `examples/common.jl`.

- The switches.
- A Maxwell-Jüttner initial state normalized with [`density`](@ref).
- The report.
- The output.

## Switches

A script takes `name=value` arguments, and every switch has a default.
The switches of `common.jl` are shared, and a script's header names the defaults it changes.

| Switches | Values | Sets |
|---|---|---|
| `grid=` | `uniform` (default), `biuniform`, `tanh`, `biuniform_270` | the momentum faces |
| `np=`, `p_max=`, `np_low=`, `p_sep=` | integers, numbers | the momentum cells and extent, and the split of the biuniform and tanh grids |
| `pitch=`, `nxi=`, `nxi_low=`, `theta_sep=` | `uniform` (default) or `theta`; integers, radians | the pitch faces |
| `T=`, `n=`, `B=` | eV, m⁻³, T | the deuterium plasma |
| `n_ar=` | m⁻³ | neutral argon added to the plasma (`relaxation.jl`) |
| `E=` or `e_over_ec=` | V/m, or multiples of [`critical_field`](@ref) at the [`thermal_coulomb_logarithm`](@ref) | the electric field |
| `mode=`, `screening=`, `lnlambda=`, `ion_size_table=` | the options of [`CollisionSettings`](@ref) | the collision frequencies |
| `synchrotron=` | `true`, `false` (default) | the synchrotron loss of `relaxation.jl` |
| `scheme=`, `p_scheme=`, `xi_scheme=` | `exponential_fitting` (default), `central`, `upwind`, `quick`, `upwind_second_order`, `van_albada`, `tcdf`, `muscl`, `smart`, `ospre` | the interpolation scheme in both directions or in one ([`AdvectionDiffusion`](@ref)) |
| `jacobian=` | `exact` (default), `frozen_limiter`, `upwind`, `central`, `exponential_fitting` | the stored Jacobian of a limiter ([`FokkerPlanckTerm`](@ref)) |
| `pmax_boundary=` | `f_zero` (default), `constant_flux`, `extrapolated_flux`, `extrapolated_flux_consistent`, `outflow` | the boundary at `p_max` ([`AdvectionDiffusion`](@ref)) |
| `nonnegative=` | `true`, `false` (default) | a [`LowerLimit`](@ref) ``f \ge 0`` on the unknown (`relaxation.jl`, `cooling.jl`) |
| `solver=` | `lu` (default), `gmres` | [`SparseLU`](@ref), or [`GMRES`](@ref) on the Newton operator |
| `preconditioner=`, `ilu_drop=` | `lu` (default), `ilu`; number | [`RadiusBlockLU`](@ref) or [`RadiusBlockILU`](@ref) for GMRES |
| `refresh=` | `every_step` (default), `every_iteration`, `on_demand` | when the preconditioner is factorized again |
| `operator=` | `exact` (default), `stored` | the Newton operator of GMRES |
| `restart=` | integer | the restart length of GMRES |
| `max_iterations=`, `residual_tolerance=`, `line_search=` | integer, number, boolean | the Newton iterations, the residual tolerance and the backtracking line search of [`NewtonSolver`](@ref) |
| `dt=`, `t_end=`, `saves=` | s, s, integer | the time points and the number of save times of [`integrate`](@ref) |
| `output=` | file name | the HDF5 output, written with [`write_output`](@ref) |
| `checkpoint=`, `checkpoint_dir=`, `checkpoint_interval=`, `overwrite=`, `resume=` | name, directory, s, boolean, boolean | the [`CheckpointSettings`](@ref) of a run of that name, and its continuation from the last checkpoint |
| `stop_file=`, `wall_budget=`, `slack=` | file name, s, s | the [`StopCondition`](@ref) of a segment |
| `log=`, `request_file=` | `stdout` or a file name; file name | the [`ProgressLog`](@ref) and the file whose appearance asks for an entry |

### Notes on the table

- `grid=biuniform_270` is a fully kinetic hot-tail grid of 270 momentum cells to ``p = 2.5``.
  It has 15 pitch cells and ignores the other grid switches.
  It is a choice of this package ([Choices of this package](choices.md)).
- `np_low=` is the number of cells below the split, and `p_sep=` its position.
- `pitch=theta` is biuniform in the pitch angle, with `nxi_low=` cells within `theta_sep=` of ``\xi = 1``.
- In `cooling.jl`, `T=` is the initial temperature.
- `pmax_boundary=extrapolated_flux` needs a `p_scheme=` other than `exponential_fitting`.
- The default scheme, exponential fitting, is a choice of this package ([Choices of this package](choices.md)).
  Its properties, the other schemes and DREAM's default are on the page [Positivity and flux limiters](model/positivity.md).

Three settings of GMRES in the examples differ from the defaults of [`GMRES`](@ref).

- The preconditioner is [`RadiusBlockLU`](@ref).
- It is factorized once per step.
- GMRES restarts after 30 iterations.

### The switches of `avalanche.jl`

`avalanche.jl` takes neither the checkpoint, stop and log switches nor the grid and plasma switches of `common.jl`.
It has switches of its own.

| Switch | Sets |
|---|---|
| `treatment=` | the time level of the knock-on term ([`avalanche.jl`](#avalanche.jl)) |
| `cutoff=` | the knock-on cutoff momentum |
| `np=`, `np_low=`, `nxi=` | the momentum and pitch cells, as in the table of [Switches](#Switches) |
| `e_over_ec=` | the field in multiples of the critical field of the total electron density |

### Lower limits

With `nonnegative=true`, a Newton step that would push an entry already at zero below it fails with `:below_lower_limit` ([`LowerLimit`](@ref), [`NewtonSolver`](@ref)).
[`FixedSteps`](@ref) then retries the step with half the time step.
At a low temperature on a wide grid, the Maxwell-Jüttner tail underflows to exact zeros, so such a run rejects and halves many steps.
Shorter steps do not remove a round-off negative update of an entry at zero.
A positive tolerance of the [`LowerLimit`](@ref) does.

### Slurm jobs

In a Slurm job with checkpoints, the segment stops before the end of the job without further switches.
[`StopCondition`](@ref) takes `SLURM_JOB_END_TIME` as its end.

## The settings of a run

[`integrate`](@ref) records the settings of the run with [`run_settings`](@ref) before the first step and returns them in `settings` of the [`TimeSolution`](@ref).
Each object of the package the run is built from has a [`describe`](@ref) record, and the record collects them by section.

- The setup: `grids`, `layout`, `system`, `updates`, `prescribed`, `solver`, `time` and `run`, each with a [`settings_hash`](@ref) and a hash of them together.
- `run` holds the `description` that a script passes to [`integrate`](@ref), such as the name of the example and the origin of its initial state.
- Outside the hashes: `code`, `environment`, `checkpoint`, `stop` and `log`.
- `code` is the [`code_version`](@ref).
  It holds the commit of Runes2P when it runs from a git checkout, as the examples do, and `"unknown"` for an installed package.
- `segment`: the number, the start, the status and the reason of a stop of this segment.

The record is taken before a stored state is restored, so restoring does not change it.
A resumed segment whose setup differs lists the changed sections in its `segment` record.
`complete` is false when some object of the run has no record of its own, such as a term defined outside the package.

## Output files

`output=run.h5` writes the solution with [`write_output`](@ref): the saved states per unknown, the layout, the solver statistics and the settings record.
[`read_output`](@ref) reads it back as a [`TimeSolution`](@ref), and [`read_settings`](@ref) reads the settings record alone.

`checkpoint=name` writes, in `checkpoint_dir=`, the files that [`CheckpointSettings`](@ref) names.

- `name.checkpoint.h5`, the last checkpoint, and `name.checkpoint.prev.h5`, the one before it.
- `name.segment<k>.h5` (`k` with four digits), the output of segment `k`.

A run that stops writes its checkpoint and exits with [`exit_status`](@ref), so that a batch script can start the next segment with `resume=true`.
[`read_output`](@ref) with the [`CheckpointSettings`](@ref) of the run joins its segments into one solution:

```julia
using Runes2P, HDF5
solution = read_output(CheckpointSettings("cooling"; directory = "run"))
f_end = unknown_view(solution.states[:, end], solution.layout, :f)   # f at the last save
```

The format of these files is on the page [Output files](model/output.md).
