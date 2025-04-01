#!/bin/bash

set -x
set -euo pipefail

rank="${OMPI_COMM_WORLD_RANK:-$SLURM_PROCID}"

[[ -z "${NSYS_FILE+x}" ]] && NSYS_FILE=report.qdrep
[[ -z "${NSYS+x}" ]] && NSYS=0

if [[ "$NSYS" -ne 0 && "$rank" -eq 0 ]]; then
  exec nsys profile \
      --cuda-um-cpu-page-faults true \
      --cuda-um-gpu-page-faults true \
      --event-sample=system-wide \
      --cpu-socket-events=61,71,265,273 \
      --cpu-socket-metrics=103,104 \
      --event-sampling-interval=1 \
      --trace=nvtx,openacc \
      --force-overwrite=true \
      -e NSYS_MPI_STORE_TEAMS_PER_RANK=1 \
      -o "$NSYS_FILE" "$@"
else
  exec "$@"
fi
