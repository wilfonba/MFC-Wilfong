# 1D air/water contact — an exact-solution gate for the semi-implicit projection method

A water slab sits in air at uniform pressure and uniform velocity in a periodic
domain, with a ~833:1 density ratio, no surface tension and no viscosity. Because
nothing drives the flow, the exact solution is that the initial profile translates
unchanged forever. Any drift in the phase densities, the pressure or the velocity is
therefore error, and its size is the metric — there is no judgement call about what
the answer should be.

```bash
./mfc.sh run examples/1D_contact_semiimplicit/case.py -n 1 -- --acfl 1 --velocity 0
./mfc.sh run examples/1D_contact_semiimplicit/case.py -n 1 -- --acfl 1 --velocity 5
./mfc.sh run examples/1D_contact_semiimplicit/case.py -n 1 -- --explicit --acfl 0.5
```

| flag | meaning |
| --- | --- |
| `--velocity 0` | stationary contact: an exact steady state, nothing may move |
| `--velocity U` | profile advects at `U` and returns to itself after `L/U` |
| `--accel A` | spatially uniform body force; accelerates both phases equally, so the exact solution survives |
| `--osc-ratio`, `--osc-freq` | oscillatory forcing `a0 + k sin(w t)`, driving the velocity through zero every half cycle |
| `--acfl` | fixed `dt` as a multiple of the water acoustic limit |
| `--explicit` | run the explicit solver as a control |

The case runs in seconds and catches in ~90 s what previously required a 35-minute
GPU run to notice.

## Acceptance criteria

Measured on the current tree, water density in pure cells (`alpha_1 > 0.99`), exact
value 1000 to within ~0.04:

| configuration | explicit | projection |
| --- | --- | --- |
| `--velocity 0` | 2.7e-09 | **0.0** |
| `--velocity 5` | 1.3e-05 | **8.7e-05** |
| `--velocity 0 --accel 100 --osc-ratio 40 --osc-freq 31` | 7.7e-04 | **3.9e+02** ← open defect |

The first two pass. The third does not, and is the open problem described below.

## The open defect

Under an oscillatory body force the water density drifts to 609–871 while the
explicit solver holds 1000.001 — a factor of 5e5. Its signature:

- **Mass is exactly conserved** (`sum(alpha_rho_k)` drifts by −0.00000000%), so nothing
  is created; `alpha_rho` and `alpha` are transported inconsistently and their ratio
  drifts.
- **The error is in the bulk, not at the interface** — it appears in cells with
  `alpha_1 > 0.99`, five or more cells away from any interface.
- **Steady acceleration is harmless** (drift 7.9e-08). The drift scales very steeply
  with the *oscillatory velocity excursion* `k/w` (12.7 m/s → 4.1e-03; 129 m/s → 390),
  so what matters is how far and how often the velocity swings.

### Mechanism

The Helmholtz right-hand side is `helm_rhs = p_adv - rho*c^2*dt*div(u*)`, and for water
`rho*c^2 = 2.64e9`. Instrumented over a 65,000-step run, the RHS non-uniformity is
essentially *all* the divergence term:

| step | `dt*div` | `rho*c^2*dt*div` | measured RHS spread |
| --- | --- | --- | --- |
| 5,400 | 1.75e-07 | 462 Pa | 632 Pa |
| 32,500 | 2.08e-02 | 5.5e+07 Pa | 4.6e+07 Pa |

**A divergence error of 1e-7 — negligible as a velocity field — becomes ~460 Pa of
spurious pressure, because it is multiplied by the bulk modulus of water.** That
pressure drives a momentum correction, which enlarges the divergence error, which is
amplified again. The loop gain is `rho*c^2`, i.e. 2.6 billion for water. This is why
the damage lives in the bulk of the stiff phase, why the explicit scheme is immune
(it has no `rho*c^2`-amplified path from a divergence error into pressure), and why
every fix aimed at the interface missed.

Note that `dt*|rhs_p_adv|` contributes only ~1.6 Pa per step at step 5,400. `rhs_p_adv`
is a *rate* (Pa/s); comparing its raw magnitude against a pressure is the mistake that
first made pressure advection look dominant.

### It is a method defect, not a porting bug

The reference implementation (`SemiImplicitFV`) was given the equivalent case
(`cases/1D_contact_bf`, with a `1D_contact_nobf` control). Its body force is
`a + b*cos(c*t + d)`, so MFC's `g + k*sin(w*t - p)` maps as `b=k, c=w, d=-pi/2`.

| | no forcing | with forcing |
| --- | --- | --- |
| reference | completes, rho 998.973–1000.117 | **NaN at step 6501** |
| MFC | rho exactly 1000.000 | survives 65,000 steps, rho 609–871 |

MFC is *better* than the reference on both, and the reference dies outright where MFC
degrades. Do not treat the reference as an oracle for this defect.

## What was tried, and why each failed

Measured against the gate. Every one of these was implemented, run, and reverted.

| attempt | result |
| --- | --- |
| **Harmonic face density** in the Helmholtz stencils | **KEPT** — see below |
| **One-sided wall correction gradient** | **KEPT** — see below |
| Rhie–Chow momentum interpolation (RHS term) | 231× better on a 1D hydrostatic column, but the 2D breakup case died at t=0.091 vs baseline 0.238 |
| Well-balanced (deviation) Rhie–Chow | died at t=0.095, same signature; the hydrostatic kink was not the killer |
| Substituting `div_u_face` in the Helmholtz RHS | broke uniform advection outright |
| Star-state upwind divergence | broke advection (ICFL 2.6e125) and produced NaNs on the gate |
| Wide Helmholtz operator matching `div∘grad` | 16× on advection, 27× on the hydrostatic interior, **no effect on the gate**; blocked by depth-1 MG ghosts |
| Energy-form solve (Avgerinos et al.), three RHS variants | all worse; exponential growth of ~2–3 per step |
| Hydrostatic wall pressure BC | made the hydrostatic column 2.7× worse |
| Mass-flux-consistent divergence via an LLF stand-in | failed ~3× earlier; needs HLLC's own face velocities |

### Two structural lessons

**`div` and `grad` compose to the WIDE Laplacian, not the compact one the solve
inverts.** The RHS divergence is a 2*dx central difference of the star velocity and
the momentum correction applies a 2*dx central pressure gradient, so their composition
is a wide operator while the solve inverts the compact 3-point one. The projection
therefore cannot enforce its own constraint. Matching the operator to `div∘grad` fixes
uniform advection by 16× and the hydrostatic interior by 27× — but does nothing for
the gate, and cannot be shipped because `mg_rho` and friends are allocated with
depth-1 ghosts (`-1:mg_m(lv)+1`) while a wide stencil needs `j±2`. Converting
multigrid would mean depth-2 ghosts everywhere, rewriting the halo pack/unpack,
undoing the depth-1 residual-slab optimization that bought −21%, and fixing RBGS,
whose `(j+k+l)`-parity colouring degenerates when `j±2` shares the centre's parity.

**Accuracy of `div(u*)` dominates its null-space properties.** The upwind face
velocity gives a compact one-sided difference with no checkerboard null space — which
was the goal — and it still blew up, because first-order truncation error is fatal
under a gain of 2.6e9. The wide central difference is very likely a deliberate
accuracy choice whose null space is the lesser evil.

### On the energy formulation

Avgerinos, Bernard, Iollo & Russo (*JCP* 393, 2019) report that solving the implicit
system **on pressure** is "more oscillatory or it diverges in the low Mach regime" on
a **collocated** mesh — an independent, published description of this exact failure —
and switch to solving on the energy. Three faithful mappings of that onto MFC were
implemented and all were worse.

The obstruction is specific to stiffened gases. For water `E ≈ 7.8e8` while
`gamma_mix*p ≈ 2.9e4`, so recovering `p` from `E` loses ~4.4 digits and seeds ~4e-7 Pa
of noise every step, which the `rho*c^2` loop then grows. The papers use gamma-law
gases where `E ~ p` and this never arises. Removing the cancellation (the property
offsets `pi_inf` and `qv` are *linear* in `alpha` and `alpha_rho`, so they cancel
analytically between `E*` and `b*`, leaving only `gamma_mix*p + ke`) improved the
static case by 2e7 — and the remaining instability still grew by 2–3 per step. The
formulation is consistent (for a stiffened gas `h/gamma_mix = c^2` exactly), so this
is a property of the mapping onto MFC's context, not of the papers' scheme.

## What is in this commit

**Harmonic face density.** The Helmholtz and multigrid stencils averaged the face
density arithmetically. Across an 833:1 jump the arithmetic mean is dominated by the
water, so the face mobility `1/rho_f` is wrong by ~3 decades on the air side; for
`div(rho^-1 grad p)` the flux-continuous coefficient is the **harmonic** mean. The
whole 1D stability map (`u` = 0…5 × ACFL 0.5…30) goes to machine precision where the
arithmetic form detonated in four steps at ACFL 30. Validated independently against
the reference, which still carries the arithmetic form and drifts by 1.14 where MFC is
now exact.

**One-sided wall correction gradient.** Next to a physical wall the pressure ghost is
a zero-order extrapolation, so the wide gradient reads a value carrying no
information — making the discrete equilibrium demand `2*rho*a` instead of `rho*a`,
which no body force can satisfy. The one-sided form is what the wide stencil returns
for a *linearly* extrapolated ghost, and needs no change to the operator's boundary
condition. On the 1D hydrostatic column: wall cells 0–9 from 120.86 to 65.64, cells
190–199 from 235.75 to 54.82. Periodic cases are bit-identical.

Both are validated on `1D_hydrostatic_semiimplicit` and this case.

## Where to go next

The projection method does not remove the low-Mach timestep restriction for
stiffened-gas multiphase, and three independent lines of evidence say the obstruction
is structural to a collocated pressure split. Before investing further:

1. **Read the limiter tag.** MFC computes five `dt` candidates (ICFL, VCFL, CCFL,
   collision, acoustic) and prints the winner. If a case is CCFL-limited, no all-Mach
   scheme helps — surface tension carries `dt ~ sqrt(rho*dx^3/sigma)`.
2. **Reduced speed of sound.** `examples/2D_interface_breakup` already softens water's
   `pi_inf` to match the gas (`water_pi_inf_mode`). Extending the same softening to
   both phases is an EOS change that fits finite-volume WENO with no structural work,
   with error `O(M^2 xi^2)`.
3. **Jacobian-free Newton–Krylov.** Solves the actual discrete system, so `div`, `grad`
   and the Laplacian never have to compose consistently — a pressure-Helmholtz
   preconditioner only has to be *approximately* right, and its inconsistency costs
   Krylov iterations rather than corrupting the answer. Everything documented above is
   fatal when the split *is* the scheme and harmless when it is only an accelerator.

Low-Mach *accuracy* is already handled: MFC carries the Thornber et al. (*JCP* 2008)
correction in HLLC and HLL via the `low_Mach` parameter.
