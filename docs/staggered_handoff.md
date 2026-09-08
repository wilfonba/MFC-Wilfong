# Staggered (MAC) discretization — handoff

Branch `jfnk`. Everything here is behind the `stagger` parameter, which defaults to `F`,
so no existing code path changes and no golden file moves.

## Why this exists

Three attempts at lifting MFC's acoustic time-step restriction have now failed on the
same root cause, and the staggered grid is the first one that addresses it directly.

1. **Semi-implicit pressure projection** (branch `projectionMethod`). Failed because a
   collocated grid cannot compose a divergence, a gradient and a Laplacian consistently:
   div and grad are both 2Δx wide, so composing them gives the *wide* Laplacian, which is
   blind to a checkerboard. A 1e-7 divergence error became hundreds of Pa through a ρc²
   gain of ~2.6e9. Full record in `examples/1D_contact_semiimplicit/README.md`.

2. **JFNK preconditioners.** Three tried (wave-speed diagonal; textbook projection
   Helmholtz; the operator consistent with what HLLC actually linearizes to). All
   measured *harmful*. The projection one fails even at ACFL 0.05 where it must reduce to
   τI. Root cause measured: `ρc²·div(z_mom/ρ)` has a gain of
   `ρ_heavy·c²_heavy/ρ_light ≈ 1.9e9` across an 833:1 interface. J has the same
   sensitivity, so M is not wrong to reproduce it — but M must then match J to about one
   part in 10⁹ to help.

3. **Staggered grid** (this document). Fixes both by construction: the operators compose
   exactly, and face velocity is a primary unknown so light-phase momentum is never
   divided by light-phase density next to a heavy cell.

## Status

| Phase | Deliverable | State |
|---|---|---|
| P0 | `stagger` param, module skeletons | done, `a9dbb634` |
| P1 | Transfers, compact div/grad, harmonic face density | **done, gate passes** |
| P2 | WENO5 to faces, conservative scalar transport | **done, order 4.90** |
| P3 | Face momentum advection | done, `2fff492b` |
| P4 | Pressure + energy, complete explicit solver | **framework done, unstable — see open problem** |
| P5 | JFNK on the staggered state | not started |
| P6 | Staggered projection preconditioner | not started |

Commits: `a9dbb634` (P0/P1), `7c6d4f8b` (P2), `2fff492b` (P3/P4), `7109bdda` (test IC fix).

Modules, all in `src/simulation/`:

- `m_staggered.fpp` — face state, cell↔face transfers, `s_face_divergence`,
  `s_cell_gradient`, `s_face_density`, `s_staggered_self_test`
- `m_weno_staggered.fpp` — WENO5-JS reconstruction to faces
- `m_rhs_staggered.fpp` — RHS assembly, time step, and the three acceptance tests

Only `num_dims == 1` is implemented.

## What is verified

All three tests run at initialization when `stagger` and `run_time_info` are on.

**P1 operator gate** — the property the whole approach turns on:

```
adjointness  |<p, div u> + <grad p, u>| / scale  =  9.77753E-17
checkerboard div(grad) min  1.60000E+05  max  1.60000E+05
checkerboard expected (compact, uniform grid)    =  1.60000E+05
```

Summation by parts holds to round-off, so div and grad are exact negative adjoints and
the pressure operator they compose is symmetric. The checkerboard response is the full
4/Δx² per dimension — where the collocated wide operator returns **exactly zero**.

**P2 scalar transport**, one full period of advection:

```
free stream  max|rhs| on a constant field   =  0.00000E+00
conservation |sum q dV - initial| / initial =  9.23706E-14
accuracy     max|q - exact| after a period  =  1.81638E-08
```

Order across resolutions (50/100/200/400 cells): errors 1.47e-5, 4.86e-7, 1.82e-8,
9.33e-10 → **4.92, 4.74, 4.28**. The fall-off at the fine end is the third-order time
integrator, not the reconstruction: halving the CFL leaves the 200-cell error nearly
unchanged (1.46e-8) while the 400-cell error drops to 4.91e-10, restoring order **4.90**.

**P3/P4 interface test** — Abgrall's condition, a two-phase 833:1 interface at uniform
pressure and velocity, one step:

```
energy-rate identity residual  1.22070E-04   (against terms of order 7.8e11)
max|p - p0| after one step     4.77012E-08
max|u - u0| after one step     2.93099E-14
```

Uniform pressure and velocity stay uniform. Note this is a statement about **one step**;
it says nothing about stability.

## The open problem

The gate diverges. At acoustic CFL 0.5 it reaches step 32 of 200. Lowering the CFL delays
but does not cure it:

| ACFL | steps reached |
|---|---|
| 0.5 | 32 |
| 0.3 | 58 |
| 0.1 | 121 |
| 0.05 | 200 (completed, but velocity had reached 41377 m/s against an exact 5.0) |

So it is instability, not CFL margin. **The scheme has no acoustic dissipation at all** —
the pressure gradient is central and the cell fluxes upwind only on the contact wave, so
the acoustic subsystem is skew-symmetric and SSP-RK3 sits on its imaginary-axis limit
(|λΔt| ≈ ACFL·π = 1.57 against √3 = 1.73 at ACFL 0.5).

### What was tried and failed

Godunov acoustic states written for the staggered layout:

```
u* = u - Δp/(2Z)     on the faces, used by the cell fluxes
p* = p - Z·Δu/2      at the cells, used by the momentum control volume
Z  = ρc from Wood's mixture speed
```

Both are provably dissipative in isolation (each contributes a positive `c·Δx/2`
Laplacian), both vanish identically for uniform p and u so every invariant above stayed
exact, and the diffusion numbers are comfortable (0.25 in water, 0.06 in air). It made
things **worse**: step 5 of 200, with either correction alone or both together.

Best current hypothesis: `u*` also becomes the momentum *advection* velocity, injecting
`∂(ρ·u·Δp/Z)/∂x` cross terms that are not dissipative. Most staggered all-speed
formulations apply the pressure-diffusion to the mass flux while treating the momentum
advection velocity separately; that distinction was not made here.

The attempt is committed as `docs/staggered_acoustic_dissipation.patch`; apply it with
`git apply docs/staggered_acoustic_dissipation.patch`. It builds and passes every
invariant — it just is not stable.

### Suggested next steps, in order

1. Work out the modified-equation form of the full nonlinear system with `u*` in place,
   rather than reasoning from the linear acoustic subsystem alone. That is where the
   non-dissipative cross terms will show up explicitly.
2. Try applying the pressure-diffusion **only to the mass and energy fluxes**, keeping the
   momentum advection velocity at `um/ρ_f`. This breaks the mass-momentum consistency that
   made velocity exact, so check the interface test immediately — it may be a real trade.
3. Failing that, an acoustic Riemann solve at the cell centres for the momentum control
   volume is the textbook route and is what the collocated path gets from HLLC.

Only after the explicit scheme is stable is P5 (JFNK on the staggered state) worth
starting — it reuses the existing Newton–GMRES verbatim, changing only pack/unpack and
the residual.

## Two design points that were hard-won

**The momentum control volume must use the mass flux the continuity equation computed.**
It is accumulated into `mfl` during the continuity transport and averaged onto the cell
centres. Reconstructing ρ independently for the momentum leaves mass and momentum
advecting on two different discrete operators, which manufactures velocity out of nothing
at a density jump — even though every equation is separately conservative and free-stream
preserving. This is what took `max|u - u0|` from 2278 m/s to 2.9e-14.

**Interface tests must set the interior and let the BCs fill the ghosts.** Writing the
profile into the ghost layers with the same predicate as the interior puts water where
periodicity requires air.

## Running it

`examples/1D_contact_staggered/case.py` is the gate: the same 833:1 air/water contact the
projection method uses, with `stagger` and `run_time_info` set.

```bash
./mfc.sh run examples/1D_contact_staggered/case.py -n 1 -- --acfl 0.05 --velocity 5
```

The three tests print at startup, before the first step, and are the reason to run this
case at all today; the run then advances with the staggered stepper and will diverge.
`--explicit` switches to the collocated solver as a control, and `--acfl` sets the step as
a multiple of the acoustic limit in water.

## Traps

Three separate times in this module a **test** failed on its own setup rather than on the
code under test:

- the operator gate's fields were not periodic in the index, so summation by parts left
  boundary terms;
- its residual was normalised by a quantity that was itself zero, because the test modes
  were orthogonal and each sum vanished independently (the operators were exact all along);
- the interface test wrote its profile into the ghost layers, so the state was
  inconsistent at the periodic seam before the solver ran.

For a scheme whose invariants are exact, verify the test's own initial condition before
believing a failure. The tell in the third case was that the worst cell was 1, not 100 —
the interface is at cell 100, so an error at cell 1 could only have come from the seam.

Two real out-of-bounds bugs were also found and fixed: the WENO face loop ran one face too
far, reading a cell past the end of the array, and the mass flux was read one face beyond
where the reconstruction can produce it.

## Related state on this branch

`jfnk` also carries two commits worth knowing about:

- `62424043` — `t_base` was never assigned, leaving `mytime` undefined for the whole JFNK
  step; `m_start_up` feeds it into a `t_stop` clamp on `dt`, so the step size was at the
  mercy of garbage. Same commit fixed the energy variable's scaling: the acoustic block
  balances only when `D_E/D_mom = γ_mix·c`, and `max|E|` was 347× too large. Krylov per
  step at ACFL 2 went 551/235/696/511/NaN → 64/99/84/101/75/66. **Any Krylov count quoted
  before this commit was measured under the corrupted clock.**
- `dfca4ce4` — a `sin(J M⁻¹v, v)` preconditioner health metric, printed under
  `run_time_info`. Judges a candidate preconditioner in one time step. Always check the
  τ→0 limit at ACFL ≈ 0.05 too, where any correct M must reduce to τI.

Outstanding and unrelated: 6 projection golden files need regenerating from the harmonic
face-density change in `0e76eef2`. They fail identically with and without any of the work
above.
