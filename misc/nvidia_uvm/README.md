## Example Workflow for NVIDIA Out-of-Core Strategy based on Unified Memory

```shell
# Allocate a node
salloc -A g183 --partition normal -t 02:00:00 -N 1 -n 4 --cpus-per-task=71

# Start uenv
uenv start --view=modules icon/25.2:v1

# cd to root directory of MFC
cd MFC-Wilfong

# Load modules
. ./mfc.sh load -c san -m g

# Build
export MFC_CUDA_CC=90
./mfc.sh build --gpu -j $(nproc) --single --unified --verbose

# Dry mfc run to compile preprocess and simulation binaries with case optimization
./mfc.sh run examples/3D_IGR_perf_test/case.py --case-optimization -t pre_process simulation --dry-run --gpu -N 1 -n 4 -j 71

# Run preprocess
./mfc.sh run examples/3D_IGR_perf_test/case.py --case-optimization -t pre_process --gpu -N 1 -n 4 -j 71

# cd to case dir
cd examples/3D_IGR_perf_test/

# run with env vars set, with binding script, and with nsys script
bash run.sh
```
The example `bind.sh`, `nsys.sh`, and `run.sh` I used can be found here under `misc/nvidia_uvm`.
These should be incorporated to the python based mfc workflow.
