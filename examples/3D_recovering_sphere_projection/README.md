# Recovering Sphere (3D, all-Mach pressure projection)

The projection counterpart of `examples/3D_recovering_sphere`: a cubic water droplet in air relaxes toward a sphere
under surface tension (sigma = 8). One octant is simulated, with symmetry planes on the low faces. The capillary flow
is slow (of order 1 m/s) against the sound speeds of water (~50 m/s here) and air (~370 m/s), so the projection steps at
the capillary limit while `--explicit` runs the HLLC solver at the acoustic one.

```shell
./mfc.sh run examples/3D_recovering_sphere_projection/case.py --case-optimization -- --N 200
```

Unlike the original, which uses the six-equation model, this case uses the five-equation model that the projection
requires, the well-balanced surface tension model and MTHINC interface compression.

## Against the explicit solver (100^3, to t = 0.07)

The projection reaches t = 0.1 in 311 steps (21 s on one A100); the explicit control needs ~37,000 steps (~18 min) to
get there and develops ~50 m/s spurious capillary currents, against the projection's 1-3 m/s, before failing at
t = 0.08. Up to t = 0.07 the interface positions agree:

| t    | axis, proj / expl | face diagonal | body diagonal |
|------|-------------------|---------------|---------------|
| 0.01 | 0.0731 / 0.0732   | 0.0976 / 0.0980 | 0.1159 / 0.1190 |
| 0.04 | 0.0745 / 0.0767   | 0.0872 / 0.0869 | 0.1000 / 0.1003 |
| 0.07 | 0.1028 / 0.1079   | 0.0913 / 0.0902 | 0.0876 / 0.0876 |

The cube's corners (body diagonal 0.130) and edges (0.106) retract while its faces (0.075) bulge past the
equivalent-sphere radius 0.093 in the first, underdamped oscillation of an inviscid droplet.
