#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
n=${1:-64}
ranks=${2:-64}
threads=${3:-1}
export OMP_NUM_THREADS=$threads MKL_NUM_THREADS=$threads
export OMP_PROC_BIND=close OMP_PLACES=cores MKL_DYNAMIC=FALSE
export OPENBLAS_NUM_THREADS=1
mpiexec -n "$ranks" ./wavelet_outer_comparison \
  -fine_n "$n" -coarse_n "$((n/8))" -level 0 -overlap 8 -ppw 20 \
  -blended_ibc_face false -mixed_direct_precision false \
  -wavelet_local_pc_precision single -compare_outers true \
  -outer_max_it 200 -outer_rtol 1e-6
