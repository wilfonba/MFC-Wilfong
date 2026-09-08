# Hysing et al. (2009) rising-bubble reference data

Finest-level results of the three benchmark groups, downloaded from
https://wwwold.mathematik.tu-dortmund.de/~featflow/en/benchmarks/cfdbenchmarking/bubble/bubble_reference.html
(`data_bench_quantities.zip`, `data_bubble_shapes.zip`) and converted to CSV.

| file prefix | group | original files |
|-------------|-------|----------------|
| `c1_tp2d`, `c2_tp2d`         | TU Dortmund (TP2D)      | `c1g1l7`, `c2g1l8` |
| `c1_freelife`, `c2_freelife` | EPFL Lausanne (FreeLIFE) | `c1g2l3`, `c2g2l3` |
| `c1_moonmd`, `c2_moonmd`     | Uni Magdeburg (MooNMD)  | `c1g3l4`, `c2g3l4` |

`c#_<group>.csv`: t, area, circularity, y_c, v_c, restricted to t <= 3 and
thinned to about 600 rows. `c#_<group>_shape.csv`: the t = 3 interface as line
segments (x0, y0, x1, y1), unthinned. Used by `../benchmark_quantities.py --case`.
