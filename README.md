# Helmholtz IBC numerical experiments

Reproduction code for the IBC manuscript's four numerical tables. This
snapshot uses continuous boundary traces, cubic normalized partitions of
unity, Gauss quadrature, and a directly evaluated residual
`norm(A*x-b)/norm(b) <= 1e-6`. Counts mean preconditioner applications.
Both right FGMRES(20) and unaccelerated iteration start from zero.

## Experiments

- Table 1: 2D, epsilon=0.01, h=k^(-2.01), delta=H/2,
  ell=ceil(0.51*log2(k)). H/h = 8,16,16,32.
- Table 2: 2D, h=k^(-3/2), delta=H/2, the same level rule.
  H/h = 8,8,8,16.
- Table 3: 2D, ppw=20, H/h=16, delta=H, ell=0,1.
- Table 4: 3D, ppw=20, H/h=8, delta=H, ell=0.

The manuscript mesh sequence is 64,128,256,512 fine elements per direction.
The archived 2D fixed-ppw log also contains historical 1024 runs; these are
not part of the current manuscript tables. No new PDE runs are needed to
inspect the saved numerical records in `results/`.

## MATLAB (2D)

MATLAB R2024a was used on Linux for the production calculations. From the
repository root, run one experiment in batch mode:

```sh
matlab -batch "addpath('matlab'); run_paper('asymptotic')"
matlab -batch "addpath('matlab'); run_paper('hthreehalves')"
matlab -batch "addpath('matlab'); run_paper('ppw20')"
```

For a small test use `run_paper('asymptotic',64)`. Dependencies are included
in `matlab/`; no Dropbox path or external project checkout is required.
Existing checkpoints are resumed. Remove or move the corresponding
`matlab/data/` output directory to start a fresh run. The 2D solves use
double precision. Drivers have a 200-application limit and 600-second
per-solve limit: a time-limited result is not equivalent to 200 iterations.

## PETSc/MUMPS (3D)

The source targets Linux with complex-double PETSc 3.25.1, MPI, Intel MKL
(including sparse BLAS headers), and MUMPS 5.9. PETSc must be built with
MUMPS and its dependencies; single-complex `cmumps` must also be available
to the compiler and linker. The production toolchain used GCC 9.5 and
MPICH 4.1.2. Set PETSC_DIR and PETSC_ARCH for your installation, with MKL
and MUMPS include/library paths available through the PETSc configuration.

```sh
make -C petsc wavelet_outer_comparison
bash petsc/run_case.sh 64 64 1
```

Run inside a scheduler allocation, not a login node. Arguments are fine
mesh size, total MPI ranks, and threads per rank. Production layouts were
64: 64x1; 128: 256x1; 256: 1024x1; 512: 1024x4. The last used 64 nodes.
The 512 case needs substantial distributed memory; do not use a workstation.
Local LU uses single precision; coarse LU, global operator, and outer
iteration use double precision. The explicit local precision option is
authoritative even if the generic console banner says double precision.

The three-dimensional level-zero space has eight continuous corner traces
per patch. The included source is the continuous-trace implementation,
not the earlier independent face-constant implementation.

## Saved results and provenance

`results/` contains numerical lines extracted from the production logs,
with machine/account metadata omitted. `source_manifest.json` records
SHA256 hashes of the original sources and logs. Solver source files are
preserved; MATLAB driver path setup is made relative for this package.
The half-overlap driver is also used for Table 1 by explicitly passing
the exponent 2.01. The manuscript itself and server credentials are not
included. Large-scale recomputation is not performed during packaging.
