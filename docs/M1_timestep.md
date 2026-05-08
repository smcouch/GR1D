# M1 Explicit Source Timestep Limiter

This note documents the runtime timestep controls for explicit M1 source
terms. The implementation is in
[`src/M1/M1_source_timestep.F90`](../src/M1/M1_source_timestep.F90), with the
final timestep selection in [`src/driver.F90`](../src/driver.F90).

The limiter applies only when M1 is active and at least one of these source
terms is explicit:

- Inelastic electron scattering, selected by `nes_evolution_type = 1`.
- Velocity-dependent energy-space coupling, selected by
  `energycoupling_evolution_type = 1`.

If both are implicit or disabled, this limiter does not affect `dt`.

## Runtime Parameters

The controls are declared in
[`src/GR1D_module.F90`](../src/GR1D_module.F90#L196-L201) and are parsed as
optional input parameters in
[`src/input_parser.F90`](../src/input_parser.F90#L189-L202):

```text
M1_source_dt_limiter
M1_source_cfl_linear
M1_source_cfl_fraction
M1_source_cfl_positive
M1_source_cfl_realizable
M1_source_dt_floor
M1_source_realizable_floor_abs
M1_source_realizable_floor_rel
M1_source_realizable_margin
M1_source_dt_cache_safety
M1_source_dt_verbose
```

Default values in the module are:

```text
M1_source_dt_limiter   = true
M1_source_cfl_linear   = 0.25d0
M1_source_cfl_fraction = 0.10d0
M1_source_cfl_positive = 0.50d0
M1_source_cfl_realizable = 0.99d0
M1_source_dt_floor     = 0.0d0
M1_source_realizable_floor_abs = 0.0d0
M1_source_realizable_floor_rel = 1.0d-12
M1_source_realizable_margin = 1.0d-8
M1_source_dt_cache_safety = 0.8d0
M1_source_dt_verbose   = 0
```

The current explicit-test parameter file uses these guardrails:

```text
M1_source_dt_limiter   = 1
M1_source_cfl_linear   = 1.0d0
M1_source_cfl_fraction = 1.0d0
M1_source_cfl_positive = 0.99d0
M1_source_cfl_realizable = 0.99d0
M1_source_dt_floor     = 1.0d-99
M1_source_realizable_floor_abs = 1.0d-12
M1_source_realizable_floor_rel = 1.0d-12
M1_source_realizable_margin = 1.0d-8
M1_source_dt_cache_safety = 0.8d0
M1_source_dt_verbose   = 1
```

## Source System Being Limited

For each active radial zone `k` and species `i`, the limiter builds a local
state vector over energy groups:

```math
U =
\left(
E_1,\ldots,E_N,
F_1,\ldots,F_N
\right)^T ,
```

where `N = number_groups`, `E_g = q_M1(k,i,g,1)`, and
`F_g = q_M1(k,i,g,2)`. This is done in
[`M1_source_timestep_limit`](../src/M1/M1_source_timestep.F90#L50-L55).

The limiter treats the explicit source update locally as

```math
\frac{dU}{dt} = S(U).
```

It evaluates separate source vectors for IES and energy coupling:

```math
S_\mathrm{ies}(U), \qquad S_\mathrm{ec}(U),
```

and, when both are enabled, also evaluates the combined source

```math
S_\mathrm{tot}(U) =
S_\mathrm{ies}(U) + S_\mathrm{ec}(U).
```

This is implemented around
[`src/M1/M1_source_timestep.F90#L57-L83`](../src/M1/M1_source_timestep.F90#L57-L83).
The reported `dt_ies` and `dt_energycoupling` are useful diagnostics, but the
actual source timestep is the minimum bound from the active combined update.

On a normal timestep, `SetTimeStep` uses source limits cached by the previous
accepted M1 explicit source update. On startup or after an invalid cache, the
direct fallback limiter updates the same supporting data used by the explicit
M1 source calculation:

```fortran
if (include_Ielectron_exp) call M1_updateeas
call M1_reconstruct
call M1_closure
```

See
[`src/M1/M1_source_timestep.F90#L46-L48`](../src/M1/M1_source_timestep.F90#L46-L48).

## Cached Source Evaluation

The expensive source algebra is normally done only in
[`M1_explicitterms`](../src/M1/M1_explicitterms.F90). The explicit IES and
energy-coupling paths compute the source rates

```math
S_\mathrm{ies}, \qquad S_\mathrm{ec}
```

and local row-sum rate bounds

```math
\lambda_\mathrm{ies}, \qquad \lambda_\mathrm{ec}.
```

Those values are reduced into a cached source timestep by
[`M1_source_timestep_cache_update`](../src/M1/M1_source_timestep.F90). The next
call to `SetTimeStep` consumes that cache instead of recomputing the source
rates. Because the cached value is one accepted step old, the timestep used by
the driver is

```math
\Delta t_\mathrm{cached}
=
\mathrm{M1\_source\_dt\_cache\_safety}
\,
\Delta t_\mathrm{source,previous}.
```

`M1_source_dt_cache_safety = 1` means no extra cache safety factor; values below
one are more conservative. If the cache is unavailable, the direct fallback in
[`M1_source_timestep_limit`](../src/M1/M1_source_timestep.F90) computes the same
kind of source rates and rate bounds immediately.

## Candidate Timestep Bounds

For each source vector and rate bound, `M1_source_dt_bound_from_lambda` computes four candidate
limits and returns their minimum:

```math
\Delta t_\mathrm{source}
= \min\left(
\Delta t_\mathrm{linear},
\Delta t_\mathrm{fraction},
\Delta t_\mathrm{positive},
\Delta t_\mathrm{realizable}
\right).
```

The implementation is in
[`M1_source_dt_bound_from_lambda`](../src/M1/M1_source_timestep.F90).

### 1. `M1_source_cfl_linear`

This controls a linearized row-sum bound. For each row of the local frozen
source operator,

```math
r_m = \sum_n |A_{mn}|,
```

and the limiter sets

```math
\lambda = \max_m r_m.
```

The candidate timestep is

```math
\Delta t_\mathrm{linear}
=
\frac{\mathrm{M1\_source\_cfl\_linear}}{\lambda}.
```

This is based on the matrix infinity norm:

```math
\rho(A) \le ||A||_\infty = \max_m \sum_n |A_{mn}|,
```

where `rho(A)` is the spectral radius. It is a conservative proxy for the
fastest linearized source timescale. It does not use the actual eigenvalues,
eigenvectors, or cancellations between matrix entries.

Interpretation:

- Smaller values make the source limiter stricter.
- Values of order `0.1` to `1` are conservative.
- Larger values, such as `10.0d0`, make this a loose guardrail rather than a
  strict stability criterion.

The code computes this at
[`src/M1/M1_source_timestep.F90#L176-L189`](../src/M1/M1_source_timestep.F90#L176-L189)
and applies the parameter at
[`src/M1/M1_source_timestep.F90#L207-L210`](../src/M1/M1_source_timestep.F90#L207-L210).

### 2. `M1_source_cfl_fraction`

This controls the maximum allowed component-wise fractional explicit source
change. The code computes

```math
\theta =
\max_m \frac{|S_m|}{q_m},
```

with scales

```math
q_m =
\max(|E_g|, f)
\quad \text{for energy rows,}
```

and

```math
q_m =
\max(|F_g|, 10^{-3}|E_g|, f)
\quad \text{for flux rows.}
```

Then

```math
\Delta t_\mathrm{fraction}
=
\frac{\mathrm{M1\_source\_cfl\_fraction}}{\theta}.
```

Equivalently, it enforces approximately

```math
\frac{|\Delta t\, S_m|}{q_m}
\le
\mathrm{M1\_source\_cfl\_fraction}
```

for every component.

This is not a strict stability theorem. It is an accuracy and nonlinear-change
guardrail. The `10^{-3}|E_g|` term in the flux scale prevents near-zero fluxes
from forcing tiny timesteps when the flux update is dynamically small compared
with the energy density.

The code computes this at
[`src/M1/M1_source_timestep.F90#L191-L199`](../src/M1/M1_source_timestep.F90#L191-L199)
and applies the parameter at
[`src/M1/M1_source_timestep.F90#L212-L214`](../src/M1/M1_source_timestep.F90#L212-L214).

### 3. `M1_source_cfl_positive`

This controls a source-only positivity estimate for the radiation energy
components. For energy rows only, if the source is decreasing the energy,

```math
S_{E_g} < 0,
```

the code computes

```math
\pi_\mathrm{pos}
=
\max_g
\left(
\frac{-S_{E_g}}{\max(E_g, f)}
\right).
```

Then

```math
\Delta t_\mathrm{positive}
=
\frac{\mathrm{M1\_source\_cfl\_positive}}{\pi_\mathrm{pos}}.
```

For a pure forward-Euler source update,

```math
E_g^{n+1} = E_g^n + \Delta t\, S_{E_g},
```

choosing `M1_source_cfl_positive <= 1` approximately prevents the source alone
from driving any energy group negative:

```math
E_g^{n+1}
\gtrsim
\left(1 - \mathrm{M1\_source\_cfl\_positive}\right) E_g^n .
```

Thus:

- `0.5` leaves roughly half the current energy under a source-only loss.
- `1.0` allows an energy group to approach zero but not cross it.
- Values above `1.0`, such as `2.0d0`, are no longer positivity preserving.
  They are only loose guardrails.

This looseness can be useful because GR1D does not apply these source pieces
as a completely isolated forward-Euler radiation update. The explicit source
contributions are subsequently handled inside the local M1 solve path, which
also has its own checks and fixes.

The code computes this at
[`src/M1/M1_source_timestep.F90#L201-L205`](../src/M1/M1_source_timestep.F90#L201-L205)
and applies the parameter at
[`src/M1/M1_source_timestep.F90#L216-L218`](../src/M1/M1_source_timestep.F90#L216-L218).

### 4. `M1_source_cfl_realizable`

This controls a CFL-like bound for the M1 realizability cone. The M1 closure
requires

```math
E_g > 0,\qquad |F_g|/X \le E_g,
```

where `X` is the radial metric factor in GR and `X = 1` in Newtonian runs.
Equivalently, both cone margins must be nonnegative:

```math
R_g^+ = E_g - F_g/X \ge 0,
```

```math
R_g^- = E_g + F_g/X \ge 0.
```

For the local source update

```math
\frac{dE_g}{dt} = S_{E_g},\qquad
\frac{dF_g}{dt} = S_{F_g},
```

the margin derivatives are

```math
\frac{dR_g^+}{dt} = S_{E_g} - S_{F_g}/X,
```

```math
\frac{dR_g^-}{dt} = S_{E_g} + S_{F_g}/X.
```

The limiter first decides whether group `g` is active enough to matter for a
global cone CFL. For one zone and species, define

```math
E_\mathrm{species} = \sum_g \max(E_g,0),
```

and

```math
E_\mathrm{active}
=
\max\left(
\mathrm{M1\_source\_realizable\_floor\_abs},
\mathrm{M1\_source\_realizable\_floor\_rel}\,
E_\mathrm{species}
\right).
```

The cone CFL is only applied to groups with

```math
E_g > E_\mathrm{active}.
```

Groups below this threshold are treated as floor-level radiation. They are not
allowed to reduce the global timestep; any tiny cone violation is left for the
local projection in `M1_implicitstep`.

For active groups, the limiter also floors the cone-margin denominator:

```math
R_\mathrm{floor}
=
\max\left(
f,\,
\mathrm{M1\_source\_realizable\_margin}\,
\max(E_g,E_\mathrm{active})
\right),
```

where `f = max(M1_source_dt_floor, tiny)`. If either margin is decreasing, the
limiter imposes

```math
\Delta t_\mathrm{realizable}
\le
\mathrm{M1\_source\_cfl\_realizable}
\frac{\max(R_g^\pm,R_\mathrm{floor})}{-dR_g^\pm/dt}.
```

This is the constraint that directly prevents the explicit IES and
energy-space advection source update from taking a realizable M1 state outside
the flux cone. Unlike the final projection in `M1_implicitstep`, this acts
before the update is taken by reducing the global timestep.

Values less than `1` retain a buffer inside the cone. The default and current
test value is

```text
M1_source_cfl_realizable = 0.99d0
```

which allows the source update to use at most 99% of the resolved cone margin.

The extra floors prevent nearly empty or already nearly free-streaming bins
from imposing a tiny global timestep just because `E_g` and `|F_g|/X` are both
very small and nearly equal.

## `M1_source_realizable_floor_abs`

This is an absolute energy-density threshold for the cone CFL. If

```math
E_g \le E_\mathrm{active},
```

and `E_active` is set by this absolute floor, the group is ignored by the
global realizability limiter. This should be used carefully because it is in
GR1D's local radiation-energy units.

The test restart currently uses

```text
M1_source_realizable_floor_abs = 1.0d-12
```

so the cone CFL does not chase floor-level radiation populations.

The same absolute floor is also used to suppress `flux>en` warning messages
from `M1_implicitstep` for floor-level bins. The projection back into the cone
still happens; only the diagnostic is suppressed below this scale.

## `M1_source_realizable_floor_rel`

This is the relative active-bin threshold. It compares each group energy to the
local species-integrated radiation energy:

```math
E_g >
\mathrm{M1\_source\_realizable\_floor\_rel}
\sum_h \max(E_h,0).
```

The default is

```text
M1_source_realizable_floor_rel = 1.0d-12
```

which only suppresses groups that are negligible compared with the local
species radiation content.

## `M1_source_realizable_margin`

This is a relative floor on the cone-margin denominator. A value of `1.0d-8`
means the global timestep limiter does not distinguish cone margins smaller
than about `10^{-8} E_g`; those tiny residual corrections are handled locally
by the post-step flux projection.

## `M1_source_dt_floor`

`M1_source_dt_floor` sets the floor `f` used in the source-timescale formulas:

```math
f =
\max(\mathrm{M1\_source\_dt\_floor}, \mathrm{tiny}).
```

It appears in the fractional-update, positivity, and realizability
denominators.

3. The positivity denominator.

4. The absolute part of the realizability margin denominator,
   `R_floor`.

A small value such as `1.0d-99` means the limiter remains sensitive to tiny
radiation populations. A larger value tells the limiter to ignore fractional
changes below that absolute scale.

## How The Source Limit Enters The Global Timestep

The normal spatial timestep is first computed in `SetTimeStep`:

```math
\Delta t_\mathrm{spatial}
=
\mathrm{dt\_reduction\_factor}\,
\mathrm{cffac}\,
\Delta t_\mathrm{hydro/M1}.
```

Then, if the source limiter is enabled,

```math
\Delta t_\mathrm{spatial}
\leftarrow
\min(\Delta t_\mathrm{spatial}, \Delta t_\mathrm{source}).
```

Finally the usual timestep-growth cap is applied:

```math
\Delta t
=
\min(\Delta t_\mathrm{spatial}, 1.05\,\Delta t_\mathrm{previous}).
```

This is implemented in
[`src/driver.F90#L72-L90`](../src/driver.F90#L72-L90).

The source limiter can only reduce the timestep. It never increases `dt`.

## Verbose Diagnostics

If

```text
M1_source_dt_verbose > 0
```

and the source limiter cuts below the spatial timestep, GR1D prints:

```text
M1 source dt limiter: dt_spatial dt_source dt_ies dt_energycoupling dt_realizable kind zone species group
```

The printed timestep values are divided by `time_gf`, matching the usual
physical-time output convention. The `kind` values are:

```text
1 = IES
2 = energy coupling
3 = combined IES + energy coupling
```

The diagnostic fields are stored for scalar output through:

```text
dt_m1_source
dt_m1_ies
dt_m1_energycoupling
dt_m1_realizable
M1_source_limiter_kind
M1_source_limiter_zone
M1_source_limiter_species
M1_source_limiter_group
```

The print statement is in
[`src/driver.F90#L84-L88`](../src/driver.F90#L84-L88).

## Practical Interpretation

The three CFL-like parameters have different meanings:

```text
M1_source_cfl_linear
```

Controls a row-sum source-operator timescale. It is closest to a linear
stability estimate, but it is conservative and not eigenvalue-sharp.

```text
M1_source_cfl_fraction
```

Controls how large the explicit source increment may be relative to the current
component scale. This is mostly an accuracy and nonlinear-change limiter.

```text
M1_source_cfl_positive
```

Controls source-only radiation energy loss. It is positivity-preserving only
for values less than or equal to about `1`.

```text
M1_source_cfl_realizable
```

Controls how much of the available M1 cone margin may be consumed by the
explicit source update. This is the direct guard against `|F|/X > E`.

For strict guardrails, use values like:

```text
M1_source_cfl_linear   = 0.25d0
M1_source_cfl_fraction = 0.10d0
M1_source_cfl_positive = 0.50d0
M1_source_cfl_realizable = 0.99d0
```

For exploratory fully explicit runs where only the row-sum and fractional
bounds are known to be too conservative, those two values can be loosened while
keeping the positivity and cone guards strict:

```text
M1_source_cfl_linear   = 10.0d0
M1_source_cfl_fraction = 5.0d0
M1_source_cfl_positive = 0.99d0
M1_source_cfl_realizable = 0.99d0
```

Those loose values should be interpreted as heuristic guardrails, not as a
mathematical proof of explicit stability.

## Limitations

The current limiter is deliberately local in radius and species. It estimates a
row-sum bound for each explicit source block, but it does not include spatial
transport, hydrodynamic feedback, or the nonlinear behavior of the full coupled
timestep.

The row-sum bound satisfies

```math
\rho(A) \le ||A||_\infty,
```

but it does not check whether the eigenvalues of `A` lie inside the actual
stability region of the time integrator. A sharper limiter would compute or
estimate the eigenvalues of the local source operator and apply the explicit
method's stability region directly, using the positivity and fractional-update
checks as additional diagnostics.
