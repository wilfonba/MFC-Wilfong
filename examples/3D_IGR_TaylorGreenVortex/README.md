# Taylor-Green Vortex (3D, IGR)

Reference:
> Hillewaert, K. (2013). TestCase C3.5 - DNS of the transition of the Taylor-Green vortex, Re=1600 - Introduction and result summary. 2nd International Workshop on high-order methods for CFD.

This is the Taylor-Green vortex case run with the information geometric regularization (IGR)
scheme, at Re = 1600.

## Sizing the problem

By default the case uses a fixed 99^3 grid. Passing `--gb-per-device` instead sizes the grid to
target a fraction (`saturation`, default 0.9) of that much memory per device, based on the number
of state scalars per cell (`scalars_per_cell`) and the working precision (`--mfc`'s `single`/`mixed`
flags, read from mfc.sh's own arguments).

## Scaling

Pass `--scaling strong` (the default) or `--scaling weak`:

- **Strong scaling**: the problem is sized once from a single device's memory budget
  (`--gb-per-device`) and stays constant as the number of devices grows.
- **Weak scaling**: the per-device problem size stays fixed and the domain grows with the number
  of devices (`nodes * tasks_per_node`, from `--mfc`). The device count is factored into `Lx x Ly x
  Lz` as close to a cube as possible, so each device still owns a perfect cube of the domain and
  the grid spacing stays uniform everywhere.

`"rdma_mpi": "F"` is the default, but if you think GPU direct MPI will work, you can try setting it to `"T"`

## Example

To run a weak-scaling case targeting ~8GB of GPU memory per device on 8 devices (2 nodes, 4 tasks
per node):

```shell
./mfc.sh run examples/3D_IGR_TaylorGreenVortex/case.py -t pre_process simulation \
             -e batch -p mypartition -N 2 -n 4 -w "01:00:00" --gpu             \
             -- --gb-per-device 8 --scaling weak
```
