# Vibrated Bubbles under a Faraday-Unstable Free Surface (2D, all-Mach pressure projection)

Two 4 mm air bubbles in a 12 cm water column under a 0.1 bar gas headspace. The closed container is vibrated vertically at
Γ = 8 and 100 Hz through MFC's oscillating body force. The free surface is far past its Faraday threshold, so its
subharmonic waves develop and interact with the breathing bubbles below. This is a proof of concept for 3D simulations.
The scaling (Bond number sets the length scale, as in `examples/2D_interface_breakup`), the forcing constraints and the
limitations are in the header of `case.py`.

```shell
./mfc.sh run examples/2D_vibrated_bubble_projection/case.py --case-optimization                   # pressure projection
./mfc.sh run examples/2D_vibrated_bubble_projection/case.py --case-optimization -- --explicit     # explicit HLLC, softened water
python3 examples/2D_vibrated_bubble_projection/analyze.py examples/2D_vibrated_bubble_projection --out result.png
```

`analyze.py` plots each bubble's depth below the free surface and its area over time, leaving out gas connected to the
headspace.

The projection keeps water's physical stiffness and steps at the flow's pace. `--explicit` must resolve the sound speed even
with water softened to the gas's 343 m/s: about 5e-8 s per step at the default 100 cells per diameter.

Limitations:

- **2D column resonance.** The bubbles are placed at h_c/2 and 2 h_c, using the 3D quasi-static critical depth for
  Bjerknes sinking, h_c = 2γp/(ρgΓ²). It does not hold in 2D. A 2D bubble is a cylinder that holds far more gas per unit
  width than a sphere, so the water above it resonates on its gas near the forcing frequency, and the pressure around the
  bubbles lags and overshoots the hydrostatic. Without the bubbles the column follows the oscillating hydrostatic pressure
  to 0.2%.
- **Volume oscillation.** The five-equation model under-predicts a bubble's volume oscillation, so the Bjerknes force comes
  out too weak.
