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
M1_source_dt_floor
M1_source_dt_verbose
```

Default values in the module are:

```text
M1_source_dt_limiter   = true
M1_source_cfl_linear   = 0.25d0
M1_source_cfl_fraction = 0.10d0
M1_source_cfl_positive = 0.50d0
M1_source_dt_floor     = 0.0d0
M1_source_dt_verbose   = 0
```

The current explicit-test parameter file uses looser guardrails:

```text
M1_source_dt_limiter   = 1
M1_source_cfl_linear   = 10.0d0
M1_source_cfl_fraction = 5.0d0
M1_source_cfl_positive = 2.0d0
M1_source_dt_floor     = 1.0d-99
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

Before evaluating these sources, the limiter updates the same supporting data
used by the explicit M1 source calculation:

```fortran
if (include_Ielectron_exp) call M1_updateeas
call M1_reconstruct
call M1_closure
```

See
[`src/M1/M1_source_timestep.F90#L46-L48`](../src/M1/M1_source_timestep.F90#L46-L48).

## Finite-Difference Jacobian

For a source vector `S(U)`, the limiter estimates the local Jacobian

```math
A_{mn} = \frac{\partial S_m}{\partial U_n}
```

by finite differences in
[`finite_difference_source`](../src/M1/M1_source_timestep.F90#L103-L138).

For component `n`, the perturbation size is

```math
h_n =
10^{-6} \max\left(|U_n|, f\right),
```

where

```math
f = \max(\mathrm{M1\_source\_dt\_floor}, \mathrm{tiny}).
```

For energy-density components close to the floor, the code uses a one-sided
difference to avoid stepping the energy through the floor:

```math
A_{mn} \approx
\frac{S_m(U + h_n e_n) - S_m(U)}{h_n}.
```

Otherwise it uses the centered difference

```math
A_{mn} \approx
\frac{S_m(U + h_n e_n) - S_m(U - h_n e_n)}{2h_n}.
```

The floor therefore affects both numerical differentiation and later
timescale denominators.

## Candidate Timestep Bounds

For each source vector and Jacobian, `source_dt_bound` computes three candidate
limits and returns their minimum:

```math
\Delta t_\mathrm{source}
= \min\left(
\Delta t_\mathrm{linear},
\Delta t_\mathrm{fraction},
\Delta t_\mathrm{positive}
\right).
```

The implementation is in
[`source_dt_bound`](../src/M1/M1_source_timestep.F90#L157-L221).

### 1. `M1_source_cfl_linear`

This controls a linearized row-sum bound. For each row of the finite-difference
Jacobian,

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

## `M1_source_dt_floor`

`M1_source_dt_floor` sets the floor `f` used in the formulas above:

```math
f =
\max(\mathrm{M1\_source\_dt\_floor}, \mathrm{tiny}).
```

It appears in three places:

1. The finite-difference perturbation scale:

   ```math
   h_n = 10^{-6}\max(|U_n|, f).
   ```

2. The fractional-update denominators.

3. The positivity denominator.

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
M1 source dt limiter: dt_spatial dt_source dt_ies dt_energycoupling kind zone species group
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

Controls a row-sum Jacobian timescale. It is closest to a linear stability
estimate, but it is conservative and not eigenvalue-sharp.

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

For strict guardrails, use values like:

```text
M1_source_cfl_linear   = 0.25d0
M1_source_cfl_fraction = 0.10d0
M1_source_cfl_positive = 0.50d0
```

For exploratory fully explicit runs where the row-sum and fractional bounds
are known to be too conservative, looser values like the current test settings
can be useful:

```text
M1_source_cfl_linear   = 10.0d0
M1_source_cfl_fraction = 5.0d0
M1_source_cfl_positive = 2.0d0
```

Those loose values should be interpreted as heuristic guardrails, not as a
mathematical proof of explicit stability.

## Limitations

The current limiter is deliberately local in radius and species. It estimates
the explicit source Jacobian for each `(zone, species)` block, but it does not
include spatial transport, hydrodynamic feedback, or the nonlinear behavior of
the full coupled timestep.

The row-sum bound satisfies

```math
\rho(A) \le ||A||_\infty,
```

but it does not check whether the eigenvalues of `A` lie inside the actual
stability region of the time integrator. A sharper limiter would compute or
estimate the eigenvalues of the local source Jacobian and apply the explicit
method's stability region directly, using the positivity and fractional-update
checks as additional diagnostics.
