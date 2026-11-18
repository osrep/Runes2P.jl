# The energy of the fully kinetic mode as the moment of f

This folder holds a prototype for issue #312, with its scripts and the logs of its runs.
The prototype takes the temperature of the cold electrons in the fully kinetic mode from the energy moment of the hot-tail distribution `f`.
It is an experiment that gives first numbers.
It is not part of the package, CI does not test it, and nobody maintains it.

## The design

In the fully kinetic mode `f` covers the bulk, and the hot region lies above the threshold `p₀`, 7 thermal momenta at `compat = :runes2p`.
The prototype replaces the row of `T_cold` of the disruption case of `test/example_disruption.jl` and keeps the rest of the case.

- The temperature: an implicit algebraic row, `0 = M_below[ε f] − 3/2 e n_cold T_cold`, where `M_below` is the [`EnergyMomentTerm`] of `f` over the region below `p₀` (`MomentumRegion` with `side = :below`).
  The envelope follows `T_cold`, so Newton's method has the column of the temperature.
  `n_cold` is the closure of the cold density, which equals the density moment of `f` below `p₀` to a relative 1e-9.
- The losses: the fluid channels of the row of the cold energy give a power density `P_fluid` per radius.
  They are the radiation and binding energy, the heat diffusion, the exchange with the ions, the heating by fluid runaways and the ohmic heating of the conductivity correction, `E²(σ − σ_kin)`.
  A cooling advection `A_p = a p (1 − w)` on `f` below the threshold applies `P_fluid`, with `w` the envelope of the hot region.
  The rate `a` per radius makes the discrete energy moment of the term equal to `P_fluid`, from the state at the start of the step.
  The term is explicit, and it moves no electrons.
- The collisions of `f` stay as they are: the Fokker-Planck term against a Maxwellian background at `T_cold`.
  The scripts measure the energy that each part of that term gives to `f`.

`E j_ohm` is not in `P_fluid`, since the field already acts on `f`.
The electron side of the ion exchange goes into an algebraic unknown that nothing reads, since the system requires both sides of the exchange.

## The runs

The case is the quench of `dream_disruption_small.h5` of the reference data, on 2 radii, with steps of 1e-7 s and the scheme `:tcdf` unless stated.
The six runs:

| Label | Run |
|---|---|
| (i) | the prototype from DREAM's state, with `T_cold` from the moment of DREAM's `f` (start (a)) |
| (ii) | the prototype from DREAM's state with its bulk replaced by a Maxwellian (start (b)) |
| (iii) | the case at `compat = :runes2p` on the row of the cold energy |
| (iv) | the case at `compat = :runes2p` on the row of the total energy |
| (v) | the case at `compat = :dream` on the row of the cold energy |
| (vi) | DREAM's values from the file |

Start (b) replaces `f` below the threshold by `c M(T*) (1 + α ξ p)`, a Maxwellian `M` shifted in pitch.
The factors keep the density of `f`, the bulk current and the moment temperature at DREAM's `T_cold`.
Run (ii') is run (ii) with the cooling that also takes back the energy of the sinks below the threshold.

## How to run the scripts

The scripts run from the root of a checkout, one Julia process at a time, with the test environment.
They need the two environment variables of the tests: `RUNES2P_ADAS_CACHE` for the ADAS cache and `RUNES2P_CHIANTI_DIR` for the CHIANTI directory.
They read the reference file through `reference_file` of `test/helpers.jl`, and write their serialized results to `output/` in this folder.

```sh
export RUNES2P_ADAS_CACHE=<ADAS cache> RUNES2P_CHIANTI_DIR=<CHIANTI directory>
D=scripts/prototypes/energy-moment
julia --project=test $D/stage1.jl                                          # starts (a) and (b), point checks
julia --project=test $D/stage2.jl 10                                       # runs (i), (ii), (ii'), (iii) to (v), 10 steps
julia --project=test $D/stage2.jl 30 pa,pb,cold                            # runs (i), (ii), (iii), 30 steps
julia --project=test $D/stage2.jl 7 total                                  # run (iv), 7 steps
PROTO_SCHEME=exponential_fitting julia --project=test $D/stage2.jl 30 pb,cold
julia --project=test $D/stage3.jl                                          # needs the three runs of stage2.jl above
julia --project=test $D/probe.jl
julia --project=test $D/probe2.jl central energy_dependent                 # also upwind, tcdf thermal, exponential_fitting
```

`stage1.jl` comes first, since the others read the starts that it writes.

| Log | Script and arguments | Commit |
|---|---|---|
| `stage1.log` | `stage1.jl` | `89056db` |
| `stage2_10.log` | `stage2.jl 10` | `89056db` |
| `stage2_30.log` | `stage2.jl 30 pa,pb,cold` | `89056db` |
| `stage2_7_total.log` | `stage2.jl 7 total` | `89056db` |
| `stage2_ef_30.log` | `stage2.jl 30 pb,cold` with `PROTO_SCHEME=exponential_fitting` | `89056db` |
| `stage3.log` | `stage3.jl` | `e51b0f0` |
| `probe.log` | `probe.jl` | `89056db` |
| `probe2.log` | `probe2.jl` with four sets of arguments, the last lines of each | `89056db` |

The scripts that wrote the logs took the path of the checkout as their first argument and wrote their results next to themselves.
The scripts here find the checkout from their own folder, write to `output/`, and read the result files under their current names.
Nothing else changed.
The lines of `probe2.log` that start with `#` were added afterwards to name the arguments of each run.

## Findings

The prototype runs: 30 steps of 1e-7 s from both starts, with 4 to 5 Newton iterations per step and no rejected step.
Its temperature lies far above that of the fluid row even from the same start.
The collisions of `f` against the fixed background give `f` energy that no other term loses, and the moment definition turns that energy into temperature.

### The temperature

`T_cold` [eV] per run, at the radii r1 / r2:

| t [s] | (i) | (ii) | (iii) | (iv) | (v) | (vi) |
|---|---|---|---|---|---|---|
| 0 | 1850 / 235.9 | 1580 / 138.8 | 1580 / 138.8 | 1580 / 138.8 | 1580 / 138.8 | 1580 / 138.8 |
| 1e-7 | 1576 / 224.5 | 1323 / 126.2 | 1316 / 121.4 | 1315 / 118.0 | 1318 / 122.9 | 1318 / 122.9 |
| 3e-7 | 1238 / 216.9 | 1015 / 114.0 | 991.9 / 101.5 | 988.8 / 93.18 | 996.5 / 102.9 | 996.5 / 102.9 |
| 7e-7 | 952.5 / 223.9 | 760.6 / 105.9 | 701.4 / 80.13 | 683.6 / 55.59 | 709.5 / 81.53 | 709.5 / 81.53 |
| 1e-6 | 869.4 / 235.7 | 685.1 / 104.8 | 596.8 / 69.81 | failed | 605.9 / 71.17 | 605.9 / 71.17 |
| 2e-6 | 793.7 / 295.1 | 604.9 / 113.3 | 420.8 / 44.79 | – | – | – |
| 3e-6 | 793.6 / 376.0 | 595.3 / 136.4 | 326.7 / 26.44 | – | – | – |

- Run (iv) fails in its eighth step, at about 7.29e-7 s, with Newton's method at its iteration limit.
- The file holds DREAM's values up to 1e-6 s.
- Start (a) gives `T_cold` 1.171 and 1.700 times DREAM's at r1 and r2.
- Run (v) differs from DREAM by at most a relative 4e-9.

The cold row from start (b) (`stage3.log`) separates the start from the evolution.
The difference is (T(ii) − T_cold row) / T_cold row:

| t [s] | r1 | r2 |
|---|---|---|
| 1e-6 | +0.149 | +0.494 |
| 3e-6 | +0.838 | +4.18 |

Under exponential fitting run (ii) ends at 519 and 65.4 eV, against 326.9 and 26.65 eV of its own cold row.

### Where the energy comes from

The energy budget of `f` closes at every step: the energy moments of its terms add up to `dW_f/dt` to 5 digits.
Under backward Euler that holds by construction, so it shows that the list of terms is complete.
The terms of run (ii), integrated over the first 10 steps [J/m³], against the excess of the bulk energy over run (iii) at 1e-6 s:

| Quantity | `:tcdf` r1 | `:tcdf` r2 | exp. fitting r1 | exp. fitting r2 |
|---|---|---|---|---|
| e-e heating of `f` | 2964 | 1132 | 2745 | 547 |
| the same on a Maxwellian at `(n_cold, T_cold)` | 923 | 246 | 153 | 2 |
| sinks | 395 | 27 | 291 | −0.5 |
| cooling (`P_fluid`) | −3870 | −402 | −3881 | −401 |
| bulk energy of (ii) − cold energy of (iii) | 3251 | 1305 | 2917 | 693 |

The e-e heating accounts for the excess energy of the bulk.
It has two parts:

- The discrete operator does not hold a Maxwellian at the temperature of its own background (table below).
- The particle source puts the new electrons into the lowest momentum cell, and the e-e collisions heat them to `T_cold`.
  The background pays nothing for that.
  In the row of the cold energy the same energy comes from the dilution of `W_cold`.

The sinks below the threshold give `f` energy too, as the row "sinks" shows.
The cooling that takes back the energy of the sinks (run (ii')) moves `T_cold` by about 1 %: 678.1 against 685.1 eV at r1 at 1e-6 s.

### The collisions on a Maxwellian

The rate at which the collisions change the energy of a Maxwellian at the background's own `(n, T)` on the grid of the case, at r1, divided by `3/2 e n T` [1/s] (`probe.log`, `probe2.log`):

| T [eV] | `:tcdf` | `:upwind` | `:central` | exponential fitting |
|---|---|---|---|---|
| 1500 | +4.7e3 | +1.8e4 | −5.1e3 | +6.2e2 |
| 500 | +4.7e4 | +8.8e4 | −6.6e4 | +9.1e3 |
| 139 | +6.7e4 | +3.2e5 | −1.8e5 | +2.3e3 |
| 50 | +1.9e5 | +1.9e6 | −3.7e5 | −1.4e1 |
| 10 | +1.8e7 | +3.3e7 | −1.8e7 | −1.1e2 |

A positive rate heats `f`.
The thermal Coulomb logarithm changes the `:tcdf` rates by less than 10 %.
With 1 % of `n` added in the lowest cell, the exponential fitting rate is +1.5e5 /s at 50 eV, the physical heating of a cold population.

### The bulk on the grid

The momentum faces lie 0.005 mₑc apart up to 0.05 mₑc and 0.0375 mₑc apart up to 0.8 mₑc.
The energy moment of a Maxwellian on this grid, divided by `3/2 n T` (`stage1.log`):

| T [eV] | 1000 | 500 | 300 | 200 | 139 | 100 | 70 | 50 to 10 | 5 |
|---|---|---|---|---|---|---|---|---|---|
| Factor | 1.010 | 1.052 | 1.031 | 0.955 | 0.951 | 0.982 | 0.998 | 1.000 | 0.944 |

- A bulk in the first momentum cell has the temperature 1.06 eV, the floor of the moment on this grid.
- The runs stay between 65 and 1850 eV, where the factor lies within 5 % of 1.
- The bulk never shrank into the first cells in these runs.

### Other observations

- The moment counts the electrons between about 3 and 7 thermal momenta as heat.
  At the end of run (ii), at r1, `p₀` lies at 0.338 mₑc.
  Below it, `f` is 101 and 3340 times the Maxwellian in the cells centred at 0.181 and 0.219 mₑc.
  A higher `T_cold` raises `p₀`, and the hot region barely forms: `n_hot` at r2 is 3.1e16 m⁻³ in run (ii) against 2.3e18 m⁻³ in run (iii) at 3e-6 s.
- Start (b) removes the seed of the hot tail, so its `n_hot`, `j_hot` and `E` are a result of the start.
  The cold row from start (b) keeps `n_hot` at r1 near 1.3e14 m⁻³ up to 1e-6 s as well.
- In run (iii) the energy moment of the bulk of `f` is 2.29 and 1.56 times `3/2 e n_cold T_cold` at 3e-6 s.
- The cooling advection alone, on a Maxwellian at r1, is shape-preserving to about 10 to 15 % in the L1 distance from the Maxwellian at its moment temperature.

| Cooling alone, 30 steps of 1e-7 s | From 1580 eV | From 139 eV |
|---|---|---|
| Rate `a` [1/s] | 1.15e5 | 5e5 |
| Final moment temperature [eV] | 1160 | 22.96 |
| L1 distance | 0.113 | 0.146 |

- Start (a) raises the conductivity with `T_cold`, so `j_ohm` jumps at the start, and `E` at r2 falls from 0.65 to 0.37 V/m in the first step of run (i).

## Shortcuts

- `P_fluid` and the envelope of the cooling are taken from the start of each step, and the cooling is explicit.
- The heat diffusion acts as a power density per radius on `f`, from the profile of `T_cold`.
- The sinks of `f` below the threshold act next to the fluid radiation.
- The ohmic heating of the conductivity correction is a fluid term, not work on `f`.
- Nothing makes the e-e collisions conserve energy.
- The particle source keeps its shape in the lowest cell.
- Start (b) is synthetic and removes the seed of the hot tail.
- The donor-cell scheme of the cooling distorts the bulk by about 10 % over a large cooling.

## What was not done

- An energy-conserving form of the self-collisions, for example a cooling that also takes back the e-e heating, was not tried.
- The adaptive stepper was not tried, since the fixed steps never failed.
- The runaway grid and the `:dream` settings of the prototype were not tried.
- No other grid was tried.
- Only `stage1.jl` was rerun after the scripts were made to run from a checkout.

[`EnergyMomentTerm`]: ../../../src/system/energy_terms.jl
