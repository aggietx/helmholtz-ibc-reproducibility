#define main hmsfem_q1_base_main
#include "hmsfem_q1_base.c"
#undef main
#include <cmumps_c.h>

/*
   Three-dimensional level-0 mIBC face multiscale coarse space and EMs-HS.

   The included HMSFEM source supplies the validated structured-Q1 Helmholtz
   element assembly and PETSc/MUMPS environment.  This file contains a separate
   coarse-space and Schwarz implementation; it does not use the polynomial
   HMSFEM coarse space.
*/

#define W3_FACES 6
#define W3_BASIS_PER_PATCH 8

typedef struct {
  PetscInt     cell, cx, cy, cz;
  PetscInt     x0, x1, y0, y1, z0, z1;
  PetscInt     nx, ny, nz, nloc;
  PetscInt    *gids;
  PetscScalar *pou;
  Mat          B;
  KSP          ksp;
  KSP          factor_ksp;
  Vec          rhs, sol;
  VecScatter   scatter;
} W3Patch;

typedef struct {
  CMUMPS_STRUC_C id;
  MUMPS_INT     *irn_loc, *jcn_loc;
  CMUMPS_COMPLEX *a_loc, *rhs_single;
  PetscInt       n, nz_loc;
  VecScatter     scatter_to_zero;
  Vec            rhs_zero;
  PetscBool      initialized, factored;
} W3CoarseMumpsSingle;

typedef struct {
  const AppCtx *ctx;
  PetscInt      coarse_n, dd_overlap, basis_solve_overlap, basis_support_overlap;
  PetscInt      block_n, nblocks, coarse_dofs;
  PetscInt      nlocal_patches;
  W3Patch      *patches;
  Mat           A, Phi, APhi, Kc;
  KSP           coarse_ksp;
  W3CoarseMumpsSingle coarse_single;
  PetscBool     use_native_single_coarse;
  Vec           coarse_rhs, coarse_sol;
  Vec           ras, residual, coarse_corr;
  PetscInt      pc_calls;
  PetscBool     blended_ibc_face;
  PetscLogDouble pc_local_solve_time, pc_coarse_solve_time;
} W3PC;

typedef struct {
  Vec       solution, residual;
  PetscReal target_abs, last_abs;
  PetscInt  max_it;
} W3TrueResidualCtx;

typedef struct {
  PetscInt     n, capacity;
  PetscInt    *rows, *cols;
  PetscScalar *vals;
} W3PBuffer;

static PetscErrorCode W3TrueResidualConverged(KSP ksp, PetscInt it, PetscReal rnorm,
                                               KSPConvergedReason *reason, void *vctx)
{
  W3TrueResidualCtx *ctx = (W3TrueResidualCtx *)vctx;
  Mat                A = NULL;
  Vec                b = NULL, solution = NULL;

  PetscFunctionBeginUser;
  (void)rnorm;
  PetscCall(KSPGetOperators(ksp, &A, NULL));
  PetscCall(KSPGetRhs(ksp, &b));
  PetscCall(KSPBuildSolution(ksp, ctx->solution, &solution));
  PetscCall(MatMult(A, solution, ctx->residual));
  PetscCall(VecAXPY(ctx->residual, -1.0, b));
  PetscCall(VecNorm(ctx->residual, NORM_2, &ctx->last_abs));
  if (PetscIsInfOrNanReal(ctx->last_abs)) *reason = KSP_DIVERGED_NANORINF;
  else if (ctx->last_abs <= ctx->target_abs) *reason = KSP_CONVERGED_RTOL;
  else if (it >= ctx->max_it) *reason = KSP_DIVERGED_ITS;
  else *reason = KSP_CONVERGED_ITERATING;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscReal W3Smoothstep(PetscReal t)
{
  t = PetscMax(0.0, PetscMin(1.0, t));
  return t * t * (3.0 - 2.0 * t);
}

static PetscReal W3RawPou1d(PetscInt node, PetscInt cell, PetscInt nc, PetscInt block_n, PetscInt overlap)
{
  PetscInt  left = cell * block_n, right = (cell + 1) * block_n;
  PetscReal wl, wr;

  if (overlap == 0) return (node >= left && node <= right) ? 1.0 : 0.0;
  wl = W3Smoothstep(((PetscReal)node - (PetscReal)(left - overlap)) / (2.0 * (PetscReal)overlap));
  wr = W3Smoothstep(((PetscReal)(right + overlap) - (PetscReal)node) / (2.0 * (PetscReal)overlap));
  if (cell == 0) return wr;
  if (cell == nc - 1) return wl;
  return PetscMin(wl, wr);
}

static PetscReal W3Pou1d(PetscInt node, PetscInt cell, PetscInt nc, PetscInt block_n, PetscInt overlap)
{
  PetscReal denom = 0.0;
  PetscInt  center = PetscMin(node / block_n, nc - 1);
  PetscInt  first = PetscMax(0, center - 2);
  PetscInt  last = PetscMin(nc - 1, center + 2);

  for (PetscInt c = first; c <= last; c++) denom += W3RawPou1d(node, c, nc, block_n, overlap);
  return denom > PETSC_SMALL ? W3RawPou1d(node, cell, nc, block_n, overlap) / denom : 0.0;
}

static inline PetscInt W3PatchNodeId(const W3Patch *p, PetscInt ix, PetscInt iy, PetscInt iz)
{
  return ((iy - p->y0) * p->nz + (iz - p->z0)) * p->nx + (ix - p->x0);
}

static void W3DecodeCell(PetscInt cell, PetscInt nc, PetscInt *cx, PetscInt *cy, PetscInt *cz)
{
  *cx = cell % nc;
  cell /= nc;
  *cz = cell % nc;
  *cy = cell / nc;
}

static PetscInt W3PatchType(PetscInt cx, PetscInt cy, PetscInt cz, PetscInt nc)
{
  PetscInt tx = cx == 0 ? 0 : (cx == nc - 1 ? 2 : 1);
  PetscInt ty = cy == 0 ? 0 : (cy == nc - 1 ? 2 : 1);
  PetscInt tz = cz == 0 ? 0 : (cz == nc - 1 ? 2 : 1);

  return tx + 3 * tz + 9 * ty;
}

static PetscErrorCode W3PBufferCreate(PetscInt capacity, W3PBuffer *buffer)
{
  PetscFunctionBeginUser;
  buffer->n = 0;
  buffer->capacity = PetscMax(capacity, 1);
  PetscCall(PetscMalloc1(buffer->capacity, &buffer->rows));
  PetscCall(PetscMalloc1(buffer->capacity, &buffer->cols));
  PetscCall(PetscMalloc1(buffer->capacity, &buffer->vals));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3PBufferAppend(W3PBuffer *buffer, PetscInt row, PetscInt col, PetscScalar value)
{
  PetscFunctionBeginUser;
  if (buffer->n == buffer->capacity) {
    PetscInt new_capacity = 2 * buffer->capacity;
    PetscCall(PetscRealloc(new_capacity * sizeof(*buffer->rows), &buffer->rows));
    PetscCall(PetscRealloc(new_capacity * sizeof(*buffer->cols), &buffer->cols));
    PetscCall(PetscRealloc(new_capacity * sizeof(*buffer->vals), &buffer->vals));
    buffer->capacity = new_capacity;
  }
  buffer->rows[buffer->n] = row;
  buffer->cols[buffer->n] = col;
  buffer->vals[buffer->n] = value;
  buffer->n++;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3PBufferDestroy(W3PBuffer *buffer)
{
  PetscFunctionBeginUser;
  PetscCall(PetscFree(buffer->vals));
  PetscCall(PetscFree(buffer->cols));
  PetscCall(PetscFree(buffer->rows));
  buffer->n = buffer->capacity = 0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3CmumpsCall(CMUMPS_STRUC_C *id)
{
  PetscFunctionBeginUser;
  PetscCall(PetscFPTrapPush(PETSC_FP_TRAP_OFF));
  cmumps_c(id);
  PetscCall(PetscFPTrapPop());
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3CoarseMumpsSingleSetup(Mat A, Vec rhs_template, W3CoarseMumpsSingle *m)
{
  PetscInt          rstart, rend, n, nz_loc = 0, k = 0;
  const PetscInt   *cols = NULL;
  const PetscScalar *vals = NULL;
  PetscInt          ncols;
  PetscMPIInt       rank;

  PetscFunctionBeginUser;
  PetscCall(PetscMemzero(m, sizeof(*m)));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCall(MatGetSize(A, &n, NULL));
  PetscCall(MatGetOwnershipRange(A, &rstart, &rend));
  for (PetscInt row = rstart; row < rend; row++) {
    PetscCall(MatGetRow(A, row, &ncols, &cols, &vals));
    nz_loc += ncols;
    PetscCall(MatRestoreRow(A, row, &ncols, &cols, &vals));
  }
  PetscCall(PetscMalloc1(nz_loc, &m->irn_loc));
  PetscCall(PetscMalloc1(nz_loc, &m->jcn_loc));
  PetscCall(PetscMalloc1(nz_loc, &m->a_loc));
  for (PetscInt row = rstart; row < rend; row++) {
    PetscCall(MatGetRow(A, row, &ncols, &cols, &vals));
    for (PetscInt j = 0; j < ncols; j++, k++) {
      m->irn_loc[k] = (MUMPS_INT)(row + 1);
      m->jcn_loc[k] = (MUMPS_INT)(cols[j] + 1);
      m->a_loc[k].r = (CMUMPS_REAL)PetscRealPart(vals[j]);
      m->a_loc[k].i = (CMUMPS_REAL)PetscImaginaryPart(vals[j]);
    }
    PetscCall(MatRestoreRow(A, row, &ncols, &cols, &vals));
  }

  m->id.par = 1;
  m->id.sym = 0;
  m->id.comm_fortran = (MUMPS_INT)MPI_Comm_c2f(PETSC_COMM_WORLD);
  m->id.job = -1;
  PetscCall(W3CmumpsCall(&m->id));
  m->initialized = PETSC_TRUE;
  PetscCheck(m->id.infog[0] >= 0, PETSC_COMM_WORLD, PETSC_ERR_LIB,
             "CMUMPS initialization failed with INFOG(1)=%d", (int)m->id.infog[0]);

  m->id.icntl[0] = -1;
  m->id.icntl[1] = -1;
  m->id.icntl[2] = 0;
  m->id.icntl[3] = 0;
  m->id.icntl[13] = 160;
  m->id.icntl[17] = 3;
  m->id.icntl[19] = 0;
  m->id.icntl[20] = 0;
  m->id.icntl[28] = 2;
  m->id.n = (MUMPS_INT)n;
  m->id.nz_loc = (MUMPS_INT)nz_loc;
  m->id.nnz_loc = (MUMPS_INT8)nz_loc;
  m->id.irn_loc = m->irn_loc;
  m->id.jcn_loc = m->jcn_loc;
  m->id.a_loc = m->a_loc;

  m->id.job = 1;
  PetscCall(W3CmumpsCall(&m->id));
  PetscCheck(m->id.infog[0] >= 0, PETSC_COMM_WORLD, PETSC_ERR_LIB,
             "CMUMPS coarse analysis failed with INFOG(1)=%d, INFOG(2)=%d",
             (int)m->id.infog[0], (int)m->id.infog[1]);
  m->id.job = 2;
  PetscCall(W3CmumpsCall(&m->id));
  PetscCheck(m->id.infog[0] >= 0, PETSC_COMM_WORLD, PETSC_ERR_LIB,
             "CMUMPS coarse factorization failed with INFOG(1)=%d, INFOG(2)=%d",
             (int)m->id.infog[0], (int)m->id.infog[1]);
  m->factored = PETSC_TRUE;

  PetscCall(VecScatterCreateToZero(rhs_template, &m->scatter_to_zero, &m->rhs_zero));
  if (!rank) PetscCall(PetscMalloc1(n, &m->rhs_single));
  m->n = n;
  m->nz_loc = nz_loc;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3CoarseMumpsSingleSolve(W3CoarseMumpsSingle *m, Vec b, Vec x)
{
  PetscMPIInt       rank;
  const PetscScalar *barray = NULL;
  PetscScalar       *xarray = NULL;

  PetscFunctionBeginUser;
  PetscCheck(m->factored, PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONGSTATE,
             "CMUMPS coarse factorization is not available");
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCall(VecScatterBegin(m->scatter_to_zero, b, m->rhs_zero, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterEnd(m->scatter_to_zero, b, m->rhs_zero, INSERT_VALUES, SCATTER_FORWARD));
  if (!rank) {
    PetscCall(VecGetArrayRead(m->rhs_zero, &barray));
    for (PetscInt i = 0; i < m->n; i++) {
      m->rhs_single[i].r = (CMUMPS_REAL)PetscRealPart(barray[i]);
      m->rhs_single[i].i = (CMUMPS_REAL)PetscImaginaryPart(barray[i]);
    }
    PetscCall(VecRestoreArrayRead(m->rhs_zero, &barray));
  }
  m->id.nrhs = 1;
  m->id.lrhs = (MUMPS_INT)m->n;
  m->id.rhs = rank ? NULL : m->rhs_single;
  m->id.job = 3;
  PetscCall(W3CmumpsCall(&m->id));
  PetscCheck(m->id.infog[0] >= 0, PETSC_COMM_WORLD, PETSC_ERR_LIB,
             "CMUMPS coarse solve failed with INFOG(1)=%d, INFOG(2)=%d",
             (int)m->id.infog[0], (int)m->id.infog[1]);
  if (!rank) {
    PetscCall(VecGetArray(m->rhs_zero, &xarray));
    for (PetscInt i = 0; i < m->n; i++) xarray[i] =
      (PetscScalar)m->rhs_single[i].r + PETSC_i * (PetscScalar)m->rhs_single[i].i;
    PetscCall(VecRestoreArray(m->rhs_zero, &xarray));
  }
  PetscCall(VecScatterBegin(m->scatter_to_zero, m->rhs_zero, x, INSERT_VALUES, SCATTER_REVERSE));
  PetscCall(VecScatterEnd(m->scatter_to_zero, m->rhs_zero, x, INSERT_VALUES, SCATTER_REVERSE));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3CoarseMumpsSingleDestroy(W3CoarseMumpsSingle *m)
{
  PetscFunctionBeginUser;
  if (m->initialized) {
    m->id.job = -2;
    PetscCall(W3CmumpsCall(&m->id));
  }
  PetscCall(VecScatterDestroy(&m->scatter_to_zero));
  PetscCall(VecDestroy(&m->rhs_zero));
  PetscCall(PetscFree(m->rhs_single));
  PetscCall(PetscFree(m->a_loc));
  PetscCall(PetscFree(m->jcn_loc));
  PetscCall(PetscFree(m->irn_loc));
  PetscCall(PetscMemzero(m, sizeof(*m)));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3SolveCoarse(W3PC *w, Vec b, Vec x)
{
  PetscFunctionBeginUser;
  if (w->use_native_single_coarse) PetscCall(W3CoarseMumpsSingleSolve(&w->coarse_single, b, x));
  else PetscCall(KSPSolve(w->coarse_ksp, b, x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscMPIInt W3OwnerFromRanges(PetscInt row, const PetscInt ranges[], PetscMPIInt size)
{
  PetscMPIInt lo = 0, hi = size;

  while (lo + 1 < hi) {
    PetscMPIInt mid = lo + (hi - lo) / 2;
    if (row < ranges[mid]) hi = mid;
    else lo = mid;
  }
  return lo;
}

static PetscErrorCode W3AssembleProlongation(Vec fine_template, PetscInt coarse_dofs, W3PBuffer *buffer, Mat *Phi)
{
  MPI_Comm          comm = PetscObjectComm((PetscObject)fine_template);
  PetscMPIInt       rank, size;
  const PetscInt   *row_ranges = NULL;
  PetscMPIInt      *send_counts = NULL, *recv_counts = NULL, *send_displs = NULL, *recv_displs = NULL, *cursor = NULL;
  PetscInt         *send_rows = NULL, *send_cols = NULL, *recv_rows = NULL, *recv_cols = NULL;
  PetscScalar      *send_vals = NULL, *recv_vals = NULL;
  PetscInt         *d_nnz = NULL, *o_nnz = NULL;
  PetscInt          rstart, rend, cstart, cend, local_rows, recv_total = 0;
  PetscMPIInt       send_total_mpi, recv_total_mpi;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(comm, &rank));
  PetscCallMPI(MPI_Comm_size(comm, &size));
  PetscCall(VecGetOwnershipRanges(fine_template, &row_ranges));
  rstart = row_ranges[rank];
  rend = row_ranges[rank + 1];
  local_rows = rend - rstart;
  cstart = (coarse_dofs / size) * rank + PetscMin((PetscInt)rank, coarse_dofs % size);
  cend = cstart + coarse_dofs / size + ((PetscInt)rank < coarse_dofs % size ? 1 : 0);

  PetscCall(PetscCalloc5(size, &send_counts, size, &recv_counts, size, &send_displs, size, &recv_displs, size, &cursor));
  for (PetscInt i = 0; i < buffer->n; i++) {
    PetscMPIInt owner = W3OwnerFromRanges(buffer->rows[i], row_ranges, size);
    send_counts[owner]++;
  }
  PetscCallMPI(MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts, 1, MPI_INT, comm));
  for (PetscMPIInt p = 1; p < size; p++) {
    send_displs[p] = send_displs[p - 1] + send_counts[p - 1];
    recv_displs[p] = recv_displs[p - 1] + recv_counts[p - 1];
  }
  if (size > 0) recv_total = recv_displs[size - 1] + recv_counts[size - 1];
  PetscCall(PetscMPIIntCast(buffer->n, &send_total_mpi));
  PetscCall(PetscMPIIntCast(recv_total, &recv_total_mpi));
  PetscCall(PetscMalloc3(PetscMax(buffer->n, 1), &send_rows, PetscMax(buffer->n, 1), &send_cols, PetscMax(buffer->n, 1), &send_vals));
  PetscCall(PetscMalloc3(PetscMax(recv_total, 1), &recv_rows, PetscMax(recv_total, 1), &recv_cols, PetscMax(recv_total, 1), &recv_vals));
  PetscCall(PetscArraycpy(cursor, send_displs, size));
  for (PetscInt i = 0; i < buffer->n; i++) {
    PetscMPIInt owner = W3OwnerFromRanges(buffer->rows[i], row_ranges, size);
    PetscMPIInt pos = cursor[owner]++;
    send_rows[pos] = buffer->rows[i];
    send_cols[pos] = buffer->cols[i];
    send_vals[pos] = buffer->vals[i];
  }
  PetscCall(W3PBufferDestroy(buffer));
  PetscCallMPI(MPI_Alltoallv(send_rows, send_counts, send_displs, MPIU_INT,
                             recv_rows, recv_counts, recv_displs, MPIU_INT, comm));
  PetscCallMPI(MPI_Alltoallv(send_cols, send_counts, send_displs, MPIU_INT,
                             recv_cols, recv_counts, recv_displs, MPIU_INT, comm));
  PetscCallMPI(MPI_Alltoallv(send_vals, send_counts, send_displs, MPIU_SCALAR,
                             recv_vals, recv_counts, recv_displs, MPIU_SCALAR, comm));
  (void)send_total_mpi;
  (void)recv_total_mpi;

  PetscCall(PetscCalloc2(PetscMax(local_rows, 1), &d_nnz, PetscMax(local_rows, 1), &o_nnz));
  for (PetscInt i = 0; i < recv_total; i++) {
    PetscInt local_row = recv_rows[i] - rstart;
    PetscCheck(local_row >= 0 && local_row < local_rows, comm, PETSC_ERR_PLIB,
               "Received prolongation row outside local ownership range");
    if (recv_cols[i] >= cstart && recv_cols[i] < cend) d_nnz[local_row]++;
    else o_nnz[local_row]++;
  }
  PetscCall(MatCreateAIJ(comm, local_rows, PETSC_DECIDE, row_ranges[size], coarse_dofs,
                         0, d_nnz, 0, o_nnz, Phi));
  PetscCall(MatSetOption(*Phi, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE));
  PetscCall(MatSetOption(*Phi, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE));
  for (PetscMPIInt src = 0; src < size; src++) {
    PetscInt begin = recv_displs[src], end = begin + recv_counts[src];
    while (begin < end) {
      PetscInt col = recv_cols[begin], next = begin + 1;
      while (next < end && recv_cols[next] == col) next++;
      PetscCall(MatSetValues(*Phi, next - begin, recv_rows + begin, 1, &col, recv_vals + begin, INSERT_VALUES));
      begin = next;
    }
  }
  PetscCall(MatAssemblyBegin(*Phi, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(*Phi, MAT_FINAL_ASSEMBLY));

  PetscCall(PetscFree2(d_nnz, o_nnz));
  PetscCall(PetscFree3(recv_rows, recv_cols, recv_vals));
  PetscCall(PetscFree3(send_rows, send_cols, send_vals));
  PetscCall(PetscFree5(send_counts, recv_counts, send_displs, recv_displs, cursor));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscReal W3FaceMassEntry(const AppCtx *ctx, PetscBool blended_face, PetscInt i, PetscInt j)
{
  PetscInt  same0 = ((i & 1) == (j & 1));
  PetscInt  same1 = (((i >> 1) & 1) == ((j >> 1) & 1));
  PetscReal mgauss = ctx->h * ctx->h * (same0 ? 2.0 : 1.0) * (same1 ? 2.0 : 1.0) / 36.0;
  PetscReal mgll = i == j ? ctx->h * ctx->h / 4.0 : 0.0;

  return blended_face ? 0.5 * (mgll + mgauss) : mgauss;
}

static PetscReal W3Level0FaceTrace(const W3Patch *p, PetscInt face, PetscInt basis,
                                   PetscInt gx, PetscInt gy, PetscInt gz)
{
  (void)face;
  PetscReal tx = (PetscReal)(gx-p->x0)/(p->x1-p->x0);
  PetscReal ty = (PetscReal)(gy-p->y0)/(p->y1-p->y0);
  PetscReal tz = (PetscReal)(gz-p->z0)/(p->z1-p->z0);
  return ((basis&1) ? tx : 1.0-tx) * ((basis&2) ? ty : 1.0-ty) *
         ((basis&4) ? tz : 1.0-tz);
}

static PetscErrorCode W3AddArtificialIbcScaled(const AppCtx *ctx, const W3Patch *p, PetscInt coarse_n, PetscBool blended_face,
                                               PetscInt ex, PetscInt ey, PetscInt ez, PetscScalar Ke[64])
{
  static const PetscInt face_nodes[6][4] = {
    {0, 2, 4, 6}, {1, 3, 5, 7},
    {0, 1, 4, 5}, {2, 3, 6, 7},
    {0, 1, 2, 3}, {4, 5, 6, 7}
  };
  PetscBool   active[6];
  PetscScalar coeff = ctx->operator_scale * ctx->abc_sign * PETSC_i /
                      (((PetscReal)ctx->ne_x / (PetscReal)coarse_n) * ctx->h);

  PetscFunctionBeginUser;
  active[0] = (PetscBool)(ex == p->x0 && p->x0 > 0);
  active[1] = (PetscBool)(ex == p->x1 - 1 && p->x1 < ctx->ne_x);
  active[2] = (PetscBool)(ey == p->y0 && p->y0 > 0);
  active[3] = (PetscBool)(ey == p->y1 - 1 && p->y1 < ctx->ne_y);
  active[4] = (PetscBool)(ez == p->z0 && p->z0 > 0);
  active[5] = (PetscBool)(ez == p->z1 - 1 && p->z1 < ctx->ne_z);
  for (PetscInt f = 0; f < 6; f++) {
      if (!active[f]) continue;
      for (PetscInt i = 0; i < 4; i++) {
        PetscInt ii = face_nodes[f][i];
        for (PetscInt j = 0; j < 4; j++) {
          PetscReal mij = W3FaceMassEntry(ctx, blended_face, i, j);
          PetscInt jj = face_nodes[f][j];
          Ke[8 * ii + jj] += coeff * mij;
        }
    }
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3BuildPatchRhs(const AppCtx *ctx, const W3Patch *p, PetscBool blended_face, Mat R)
{
  static const PetscInt face_nodes[6][4] = {
    {0, 2, 4, 6}, {1, 3, 5, 7},
    {0, 1, 4, 5}, {2, 3, 6, 7},
    {0, 1, 2, 3}, {4, 5, 6, 7}
  };

  PetscFunctionBeginUser;
  PetscCall(MatZeroEntries(R));
  for (PetscInt ey = p->y0; ey < p->y1; ey++) {
    for (PetscInt ez = p->z0; ez < p->z1; ez++) {
      for (PetscInt ex = p->x0; ex < p->x1; ex++) {
        PetscBool active[6] = {
          (PetscBool)(ex == p->x0), (PetscBool)(ex == p->x1 - 1),
          (PetscBool)(ey == p->y0), (PetscBool)(ey == p->y1 - 1),
          (PetscBool)(ez == p->z0), (PetscBool)(ez == p->z1 - 1)
        };
        PetscInt rows[8] = {
          W3PatchNodeId(p, ex,     ey,     ez),
          W3PatchNodeId(p, ex + 1, ey,     ez),
          W3PatchNodeId(p, ex,     ey + 1, ez),
          W3PatchNodeId(p, ex + 1, ey + 1, ez),
          W3PatchNodeId(p, ex,     ey,     ez + 1),
          W3PatchNodeId(p, ex + 1, ey,     ez + 1),
          W3PatchNodeId(p, ex,     ey + 1, ez + 1),
          W3PatchNodeId(p, ex + 1, ey + 1, ez + 1)
        };
        PetscInt gx[8] = {ex, ex + 1, ex, ex + 1, ex, ex + 1, ex, ex + 1};
        PetscInt gy[8] = {ey, ey, ey + 1, ey + 1, ey, ey, ey + 1, ey + 1};
        PetscInt gz[8] = {ez, ez, ez, ez, ez + 1, ez + 1, ez + 1, ez + 1};
        for (PetscInt f = 0; f < W3_FACES; f++) {
          if (!active[f]) continue;
          for (PetscInt q = 0; q < W3_BASIS_PER_PATCH; q++) {
            PetscInt col = q;
            for (PetscInt i = 0; i < 4; i++) {
              PetscScalar v = 0.0;
              for (PetscInt j = 0; j < 4; j++) {
                PetscInt node = face_nodes[f][j];
                v += W3FaceMassEntry(ctx, blended_face, i, j) *
                     W3Level0FaceTrace(p, f, q, gx[node], gy[node], gz[node]);
              }
              PetscCall(MatSetValue(R, rows[face_nodes[f][i]], col, v, ADD_VALUES));
            }
          }
        }
      }
    }
  }
  PetscCall(MatAssemblyBegin(R, MAT_FINAL_ASSEMBLY));
  PetscCall(MatAssemblyEnd(R, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3PatchSetup(const AppCtx *ctx, W3PC *w, W3Patch *p, PetscInt cell,
                                   PetscInt patch_overlap, PetscInt pou_overlap,
                                   Vec global_template, KSP shared_ksp, PetscBool build_runtime, PetscBool build_basis,
                                   W3PBuffer *pbuffer, PetscReal *basis_min_norm,
                                   PetscLogDouble *local_solve_time, PetscLogDouble *p_buffer_time)
{
  PetscInt       rows[8], *to = NULL;
  Mat            R = NULL, X = NULL, F = NULL;
  const PetscScalar *xarr = NULL;
  IS             is_from = NULL, is_to = NULL;
  PC             pc;
  PetscLogDouble t0, t1, t2;

  PetscFunctionBeginUser;
  PetscCall(PetscTime(&t0));
  PetscCall(PetscMemzero(p, sizeof(*p)));
  p->cell = cell;
  W3DecodeCell(cell, w->coarse_n, &p->cx, &p->cy, &p->cz);
  p->x0 = PetscMax(0, p->cx * w->block_n - patch_overlap);
  p->x1 = PetscMin(ctx->ne_x, (p->cx + 1) * w->block_n + patch_overlap);
  p->y0 = PetscMax(0, p->cy * w->block_n - patch_overlap);
  p->y1 = PetscMin(ctx->ne_y, (p->cy + 1) * w->block_n + patch_overlap);
  p->z0 = PetscMax(0, p->cz * w->block_n - patch_overlap);
  p->z1 = PetscMin(ctx->ne_z, (p->cz + 1) * w->block_n + patch_overlap);
  p->nx = p->x1 - p->x0 + 1;
  p->ny = p->y1 - p->y0 + 1;
  p->nz = p->z1 - p->z0 + 1;
  p->nloc = p->nx * p->ny * p->nz;

  PetscCall(PetscMalloc1(p->nloc, &p->gids));
  PetscCall(PetscMalloc1(p->nloc, &p->pou));
  PetscCall(PetscMalloc1(p->nloc, &to));
  for (PetscInt iy = p->y0; iy <= p->y1; iy++) {
    for (PetscInt iz = p->z0; iz <= p->z1; iz++) {
      for (PetscInt ix = p->x0; ix <= p->x1; ix++) {
        PetscInt lid = W3PatchNodeId(p, ix, iy, iz);
        p->gids[lid] = NodeId(ctx, ix, iy, iz);
        p->pou[lid] = W3Pou1d(ix, p->cx, w->coarse_n, w->block_n, pou_overlap) *
                      W3Pou1d(iy, p->cy, w->coarse_n, w->block_n, pou_overlap) *
                      W3Pou1d(iz, p->cz, w->coarse_n, w->block_n, pou_overlap);
        to[lid] = lid;
      }
    }
  }

  if (!shared_ksp) {
    PetscCall(MatCreateSeqAIJ(PETSC_COMM_SELF, p->nloc, p->nloc, 27, NULL, &p->B));
    PetscCall(MatSetOption(p->B, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
    for (PetscInt ey = p->y0; ey < p->y1; ey++) {
      for (PetscInt ez = p->z0; ez < p->z1; ez++) {
        for (PetscInt ex = p->x0; ex < p->x1; ex++) {
          PetscScalar Ke[64];
          PetscCall(ComputeElementMatrixOnly(ctx, ex, ey, ez, Ke));
          PetscCall(W3AddArtificialIbcScaled(ctx, p, w->coarse_n, w->blended_ibc_face, ex, ey, ez, Ke));
          rows[0] = W3PatchNodeId(p, ex,     ey,     ez);
          rows[1] = W3PatchNodeId(p, ex + 1, ey,     ez);
          rows[2] = W3PatchNodeId(p, ex,     ey + 1, ez);
          rows[3] = W3PatchNodeId(p, ex + 1, ey + 1, ez);
          rows[4] = W3PatchNodeId(p, ex,     ey,     ez + 1);
          rows[5] = W3PatchNodeId(p, ex + 1, ey,     ez + 1);
          rows[6] = W3PatchNodeId(p, ex,     ey + 1, ez + 1);
          rows[7] = W3PatchNodeId(p, ex + 1, ey + 1, ez + 1);
          PetscCall(MatSetValues(p->B, 8, rows, 8, rows, Ke, ADD_VALUES));
        }
      }
    }
    PetscCall(MatAssemblyBegin(p->B, MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(p->B, MAT_FINAL_ASSEMBLY));
  }

  if (build_runtime) {
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, p->nloc, &p->rhs));
    PetscCall(VecDuplicate(p->rhs, &p->sol));
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, p->nloc, p->gids, PETSC_COPY_VALUES, &is_from));
    PetscCall(ISCreateGeneral(PETSC_COMM_SELF, p->nloc, to, PETSC_OWN_POINTER, &is_to));
    to = NULL;
    PetscCall(VecScatterCreate(global_template, is_from, p->rhs, is_to, &p->scatter));
    PetscCall(ISDestroy(&is_from));
    PetscCall(ISDestroy(&is_to));
  }
  PetscCall(PetscFree(to));

  if (!shared_ksp) {
    PetscCall(KSPCreate(PETSC_COMM_SELF, &p->ksp));
    PetscCall(KSPSetOperators(p->ksp, p->B, p->B));
    PetscCall(KSPSetType(p->ksp, KSPPREONLY));
    PetscCall(KSPGetPC(p->ksp, &pc));
    PetscCall(PCSetType(pc, PCLU));
    PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    PetscCall(KSPSetOptionsPrefix(p->ksp, "wavelet_local_"));
    PetscCall(KSPSetFromOptions(p->ksp));
    PetscCall(KSPSetUp(p->ksp));
    p->factor_ksp = p->ksp;
  } else {
    Mat      shared_A = NULL;
    PetscInt shared_m, shared_n;

    PetscCall(KSPGetOperators(shared_ksp, &shared_A, NULL));
    PetscCall(MatGetSize(shared_A, &shared_m, &shared_n));
    PetscCheck(shared_m == p->nloc && shared_n == p->nloc, PETSC_COMM_SELF, PETSC_ERR_ARG_SIZ,
               "Shared patch factor has incompatible dimensions");
    p->factor_ksp = shared_ksp;
  }

  if (build_basis) {
    PetscCall(MatCreateSeqDense(PETSC_COMM_SELF, p->nloc, W3_BASIS_PER_PATCH, NULL, &R));
    PetscCall(MatDuplicate(R, MAT_DO_NOT_COPY_VALUES, &X));
    PetscCall(W3BuildPatchRhs(ctx, p, w->blended_ibc_face, R));
    PetscCall(KSPGetPC(p->factor_ksp, &pc));
    PetscCall(PCFactorGetMatrix(pc, &F));
    PetscCall(MatMatSolve(F, R, X));
    PetscCall(MatDenseGetArrayRead(X, &xarr));
    PetscCall(PetscTime(&t1));
    for (PetscInt j = 0; j < W3_BASIS_PER_PATCH; j++) {
      PetscInt  gcol = p->cell * W3_BASIS_PER_PATCH + j;
      PetscReal n2 = 0.0;
      for (PetscInt i = 0; i < p->nloc; i++) {
        PetscScalar value = p->pou[i] * xarr[(size_t)j * (size_t)p->nloc + (size_t)i];
        n2 += PetscAbsScalar(value) * PetscAbsScalar(value);
        if (p->pou[i] != 0.0) PetscCall(W3PBufferAppend(pbuffer, p->gids[i], gcol, value));
      }
      *basis_min_norm = PetscMin(*basis_min_norm, PetscSqrtReal(n2));
    }
    PetscCall(MatDenseRestoreArrayRead(X, &xarr));
    PetscCall(PetscTime(&t2));
    *local_solve_time += t1 - t0;
    *p_buffer_time += t2 - t1;
    PetscCall(MatDestroy(&X));
    PetscCall(MatDestroy(&R));
  } else {
    PetscCall(PetscTime(&t1));
    *local_solve_time += t1 - t0;
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3PatchDestroy(W3Patch *p)
{
  PetscFunctionBeginUser;
  PetscCall(VecScatterDestroy(&p->scatter));
  PetscCall(VecDestroy(&p->sol));
  PetscCall(VecDestroy(&p->rhs));
  PetscCall(KSPDestroy(&p->ksp));
  PetscCall(MatDestroy(&p->B));
  PetscCall(PetscFree(p->pou));
  PetscCall(PetscFree(p->gids));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3CheckPou(const AppCtx *ctx, PetscInt coarse_n, PetscInt block_n, PetscInt overlap, PetscReal *defect)
{
  PetscReal local = 0.0, global = 0.0;
  PetscMPIInt rank, size;
  PetscInt ndof = (ctx->ne_x + 1) * (ctx->ne_y + 1) * (ctx->ne_z + 1);

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  for (PetscInt gid = rank; gid < ndof; gid += size) {
    PetscInt ix = gid % (ctx->ne_x + 1);
    PetscInt q = gid / (ctx->ne_x + 1);
    PetscInt iz = q % (ctx->ne_z + 1);
    PetscInt iy = q / (ctx->ne_z + 1);
    PetscReal sx = 0.0, sy = 0.0, sz = 0.0;

    for (PetscInt c = 0; c < coarse_n; c++) {
      sx += W3Pou1d(ix, c, coarse_n, block_n, overlap);
      sy += W3Pou1d(iy, c, coarse_n, block_n, overlap);
      sz += W3Pou1d(iz, c, coarse_n, block_n, overlap);
    }
    local = PetscMax(local, PetscAbsReal(sx * sy * sz - 1.0));
  }
  PetscCallMPI(MPI_Allreduce(&local, &global, 1, MPIU_REAL, MPI_MAX, PETSC_COMM_WORLD));
  *defect = global;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3Setup(const AppCtx *ctx, Mat A, Vec b, PetscInt coarse_n, PetscInt dd_overlap,
                              PetscInt basis_solve_overlap, PetscInt basis_support_overlap,
                              PetscBool blended_ibc_face,
                              PetscBool use_native_single_coarse, W3PC *w,
                              PetscLogDouble *local_setup, PetscLogDouble *galerkin_setup, PetscLogDouble *coarse_factor)
{
  PetscMPIInt    rank, size;
  PetscReal      local_min = PETSC_MAX_REAL, global_min = 0.0, dd_pou_defect = 0.0, basis_pou_defect = 0.0;
  PetscLogDouble t0, t_insert, t1, t_ap, t_ht, t2, t3;
  PetscLogDouble dd_setup_time = 0.0, basis_setup_time = 0.0, p_buffer_time = 0.0;
  PetscLogDouble dd_setup_max = 0.0, basis_setup_max = 0.0, p_buffer_max = 0.0;
  Mat            PhiH = NULL;
  PC             pc;
  W3PBuffer      pbuffer;
  W3Patch       *basis_patches = NULL;
  PetscBool      full_shared_patch_set = (PetscBool)(dd_overlap == basis_solve_overlap &&
                                                      dd_overlap == basis_support_overlap);
  PetscBool      reuse_dd_factors = (PetscBool)(!full_shared_patch_set &&
                                                 dd_overlap == basis_solve_overlap);
  PetscInt       patch_width, max_entries_per_patch;
  PetscInt       dd_factor_owner[27], basis_factor_owner[27];
  PetscInt       dd_factor_count = 0, basis_factor_count = 0;
  PetscInt       dd_factor_min = 0, dd_factor_max = 0, basis_factor_min = 0, basis_factor_max = 0;

  PetscFunctionBeginUser;
  PetscCall(PetscMemzero(w, sizeof(*w)));
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  w->ctx = ctx;
  w->A = A;
  w->coarse_n = coarse_n;
  w->dd_overlap = dd_overlap;
  w->basis_solve_overlap = basis_solve_overlap;
  w->basis_support_overlap = basis_support_overlap;
  w->blended_ibc_face = blended_ibc_face;
  w->use_native_single_coarse = use_native_single_coarse;
  w->block_n = ctx->ne_x / coarse_n;
  w->nblocks = coarse_n * coarse_n * coarse_n;
  w->coarse_dofs = W3_BASIS_PER_PATCH * w->nblocks;
  PetscCheck(w->nblocks % (PetscInt)size == 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
             "Current VecScatter schedule requires coarse_n^3 divisible by MPI size");
  w->nlocal_patches = w->nblocks / (PetscInt)size;
  PetscCall(PetscCalloc1(w->nlocal_patches, &w->patches));
  if (!full_shared_patch_set) PetscCall(PetscCalloc1(w->nlocal_patches, &basis_patches));
  PetscCall(W3CheckPou(ctx, coarse_n, w->block_n, dd_overlap, &dd_pou_defect));
  PetscCall(W3CheckPou(ctx, coarse_n, w->block_n, basis_support_overlap, &basis_pou_defect));
  PetscCheck(dd_pou_defect < 1e-6, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
             "DD cubic POU defect %g", (double)dd_pou_defect);
  PetscCheck(basis_pou_defect < 1e-6, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
             "Basis cubic POU defect %g", (double)basis_pou_defect);

  patch_width = PetscMin(ctx->ne_x, w->block_n + 2 * basis_support_overlap) + 1;
  max_entries_per_patch = patch_width * patch_width * patch_width * W3_BASIS_PER_PATCH;
  PetscCall(W3PBufferCreate(w->nlocal_patches * max_entries_per_patch, &pbuffer));
  for (PetscInt type = 0; type < 27; type++) {
    dd_factor_owner[type] = -1;
    basis_factor_owner[type] = -1;
  }
  PetscCall(PetscTime(&t0));

  if (!full_shared_patch_set) {
    for (PetscInt j = 0; j < w->nlocal_patches; j++) {
      PetscInt cell = (PetscInt)rank + j * (PetscInt)size;
      PetscInt cx, cy, cz, type;
      KSP      shared_ksp = NULL;

      W3DecodeCell(cell, w->coarse_n, &cx, &cy, &cz);
      type = W3PatchType(cx, cy, cz, w->coarse_n);
      if (dd_factor_owner[type] >= 0) shared_ksp = w->patches[dd_factor_owner[type]].factor_ksp;
      PetscCall(W3PatchSetup(ctx, w, &w->patches[j], cell, dd_overlap, dd_overlap, b, shared_ksp,
                             PETSC_TRUE, PETSC_FALSE, NULL, NULL, &dd_setup_time, NULL));
      if (dd_factor_owner[type] < 0) {
        dd_factor_owner[type] = j;
        dd_factor_count++;
      }
    }
  }
  for (PetscInt j = 0; j < w->nlocal_patches; j++) {
    PetscInt cell = (PetscInt)rank + j * (PetscInt)size;
    PetscInt cx, cy, cz, type;
    KSP      shared_ksp = NULL;
    W3Patch *p = full_shared_patch_set ? &w->patches[j] : &basis_patches[j];

    W3DecodeCell(cell, w->coarse_n, &cx, &cy, &cz);
    type = W3PatchType(cx, cy, cz, w->coarse_n);
    if (reuse_dd_factors) {
      shared_ksp = w->patches[j].factor_ksp;
    } else if (basis_factor_owner[type] >= 0) {
      W3Patch *owner = full_shared_patch_set ? &w->patches[basis_factor_owner[type]] : &basis_patches[basis_factor_owner[type]];
      shared_ksp = owner->factor_ksp;
    }
    PetscCall(W3PatchSetup(ctx, w, p, cell, basis_solve_overlap, basis_support_overlap, b, shared_ksp,
                           full_shared_patch_set, PETSC_TRUE, &pbuffer, &local_min,
                           &basis_setup_time, &p_buffer_time));
    if (!reuse_dd_factors && basis_factor_owner[type] < 0) {
      basis_factor_owner[type] = j;
      basis_factor_count++;
    }
    if (reuse_dd_factors) PetscCall(W3PatchDestroy(p));
  }
  PetscCall(PetscTime(&t_insert));
  PetscCall(W3AssembleProlongation(b, w->coarse_dofs, &pbuffer, &w->Phi));
  PetscCall(PetscTime(&t1));
  if (!full_shared_patch_set) {
    if (!reuse_dd_factors) {
      for (PetscInt j = 0; j < w->nlocal_patches; j++) PetscCall(W3PatchDestroy(&basis_patches[j]));
    }
    PetscCall(PetscFree(basis_patches));
  }
  PetscCallMPI(MPI_Reduce(&dd_setup_time, &dd_setup_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&basis_setup_time, &basis_setup_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&p_buffer_time, &p_buffer_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&dd_factor_count, &dd_factor_min, 1, MPIU_INT, MPI_MIN, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&dd_factor_count, &dd_factor_max, 1, MPIU_INT, MPI_MAX, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&basis_factor_count, &basis_factor_min, 1, MPIU_INT, MPI_MIN, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&basis_factor_count, &basis_factor_max, 1, MPIU_INT, MPI_MAX, 0, PETSC_COMM_WORLD));
  if (full_shared_patch_set) {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "setup phase shared DD/basis matrices, LU, and six-RHS solves (max rank): %.3f s\n", (double)basis_setup_max));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "shared DD/basis LU factors per rank (min/max): %" PetscInt_FMT "/%" PetscInt_FMT "\n",
                          basis_factor_min, basis_factor_max));
  } else {
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "setup phase DD local matrices and LU (max rank): %.3f s\n", (double)dd_setup_max));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "DD local LU factors per rank (min/max): %" PetscInt_FMT "/%" PetscInt_FMT "\n",
                          dd_factor_min, dd_factor_max));
    if (reuse_dd_factors) {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "setup phase basis six-RHS solves and restriction using DD LU (max rank): %.3f s\n", (double)basis_setup_max));
    } else {
      PetscCall(PetscPrintf(PETSC_COMM_WORLD,
        "setup phase basis local matrices, LU, and six-RHS solves (max rank): %.3f s\n", (double)basis_setup_max));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "basis local LU factors per rank (min/max): %" PetscInt_FMT "/%" PetscInt_FMT "\n",
                            basis_factor_min, basis_factor_max));
    }
  }
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "setup phase P buffer fill (max rank): %.3f s\n", (double)p_buffer_max));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "setup phase distributed P assembly: %.3f s\n", (double)(t1 - t_insert)));
  PetscCallMPI(MPI_Allreduce(&local_min, &global_min, 1, MPIU_REAL, MPI_MIN, PETSC_COMM_WORLD));
  PetscCheck(global_min > 100.0 * PETSC_MACHINE_EPSILON, PETSC_COMM_WORLD, PETSC_ERR_PLIB,
             "Wavelet coarse basis contains a zero or numerically empty column; minimum column norm %g", (double)global_min);

  PetscCall(MatMatMult(A, w->Phi, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &w->APhi));
  PetscCall(PetscTime(&t_ap));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "setup phase AP: %.3f s\n", (double)(t_ap - t1)));
  PetscCall(MatHermitianTranspose(w->Phi, MAT_INITIAL_MATRIX, &PhiH));
  PetscCall(PetscTime(&t_ht));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "setup phase P Hermitian transpose: %.3f s\n", (double)(t_ht - t_ap)));
  PetscCall(MatMatMult(PhiH, w->APhi, MAT_INITIAL_MATRIX, PETSC_DEFAULT, &w->Kc));
  PetscCall(MatDestroy(&PhiH));
  PetscCall(MatDestroy(&w->APhi));
  PetscCall(PetscTime(&t2));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "setup phase P^*AP: %.3f s\n", (double)(t2 - t_ht)));
  PetscCall(MatCreateVecs(w->Kc, &w->coarse_sol, &w->coarse_rhs));
  PetscCall(VecDuplicate(b, &w->ras));
  PetscCall(VecDuplicate(b, &w->residual));
  PetscCall(VecDuplicate(b, &w->coarse_corr));

  if (w->use_native_single_coarse) {
    PetscCall(W3CoarseMumpsSingleSetup(w->Kc, w->coarse_rhs, &w->coarse_single));
  } else {
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &w->coarse_ksp));
    PetscCall(KSPSetOperators(w->coarse_ksp, w->Kc, w->Kc));
    PetscCall(KSPSetType(w->coarse_ksp, KSPPREONLY));
    PetscCall(KSPGetPC(w->coarse_ksp, &pc));
    PetscCall(PCSetType(pc, PCLU));
    PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
    PetscCall(KSPSetOptionsPrefix(w->coarse_ksp, "wavelet_coarse_"));
    PetscCall(KSPSetFromOptions(w->coarse_ksp));
    PetscCall(KSPSetUp(w->coarse_ksp));
  }
  PetscCall(PetscTime(&t3));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "setup phase coarse LU: %.3f s\n", (double)(t3 - t2)));
  *local_setup = t1 - t0;
  *galerkin_setup = t2 - t1;
  *coarse_factor = t3 - t2;
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "level-0 space: continuous corner traces lifted with modified IBC, basis/patch=8, coarse dofs=%" PetscInt_FMT ", basis solve/support delta/h=%" PetscInt_FMT "/%" PetscInt_FMT ", basis cubic POU defect=%.3e, DD cubic POU defect=%.3e, minimum basis column norm=%.3e\n",
    w->coarse_dofs, basis_solve_overlap, basis_support_overlap,
    (double)basis_pou_defect, (double)dd_pou_defect, (double)global_min));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3ApplyRas(W3PC *w, Vec b, Vec x)
{
  const PetscScalar *u = NULL;

  PetscFunctionBeginUser;
  PetscCall(VecSet(x, 0.0));
  for (PetscInt j = 0; j < w->nlocal_patches; j++) {
    W3Patch *p = &w->patches[j];
    PetscCall(VecScatterBegin(p->scatter, b, p->rhs, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(VecScatterEnd(p->scatter, b, p->rhs, INSERT_VALUES, SCATTER_FORWARD));
    PetscCall(KSPSolve(p->factor_ksp, p->rhs, p->sol));
    PetscCall(VecGetArrayRead(p->sol, &u));
    for (PetscInt i = 0; i < p->nloc; i++) {
      if (p->pou[i] != 0.0) PetscCall(VecSetValue(x, p->gids[i], p->pou[i] * u[i], ADD_VALUES));
    }
    PetscCall(VecRestoreArrayRead(p->sol, &u));
  }
  PetscCall(VecAssemblyBegin(x));
  PetscCall(VecAssemblyEnd(x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3PCApply(PC pc_shell, Vec b, Vec x)
{
  W3PC          *w = NULL;
  void          *ptr = NULL;
  PetscLogDouble t0, t1, t2;

  PetscFunctionBeginUser;
  PetscCall(PCShellGetContext(pc_shell, &ptr));
  w = (W3PC *)ptr;
  PetscCall(PetscTime(&t0));
  PetscCall(W3ApplyRas(w, b, w->ras));
  PetscCall(PetscTime(&t1));
  PetscCall(MatMult(w->A, w->ras, w->residual));
  PetscCall(VecAYPX(w->residual, -1.0, b)); /* residual = b - A B_RAS b */
  PetscCall(MatMultHermitianTranspose(w->Phi, w->residual, w->coarse_rhs));
  PetscCall(W3SolveCoarse(w, w->coarse_rhs, w->coarse_sol));
  PetscCall(MatMult(w->Phi, w->coarse_sol, w->coarse_corr));
  PetscCall(VecWAXPY(x, 1.0, w->coarse_corr, w->ras));
  PetscCall(PetscTime(&t2));
  w->pc_local_solve_time += t1 - t0;
  w->pc_coarse_solve_time += t2 - t1;
  w->pc_calls++;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3CoarseSolution(W3PC *w, Vec b, Vec x, PetscReal *relres, PetscLogDouble *solve_time)
{
  PetscReal      bn, rn;
  PetscLogDouble t0, t1;

  PetscFunctionBeginUser;
  PetscCall(PetscTime(&t0));
  PetscCall(MatMultHermitianTranspose(w->Phi, b, w->coarse_rhs));
  PetscCall(W3SolveCoarse(w, w->coarse_rhs, w->coarse_sol));
  PetscCall(MatMult(w->Phi, w->coarse_sol, x));
  PetscCall(PetscTime(&t1));
  PetscCall(MatMult(w->A, x, w->residual));
  PetscCall(VecAXPY(w->residual, -1.0, b));
  PetscCall(VecNorm(w->residual, NORM_2, &rn));
  PetscCall(VecNorm(b, NORM_2, &bn));
  *relres = rn / bn;
  *solve_time = t1 - t0;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3RelativeFeErrors(const AppCtx *ctx, Vec u, Vec ref, PetscReal *l2rel, PetscReal *hk1rel)
{
  VecScatter      su = NULL, sr = NULL;
  Vec             ua = NULL, ra = NULL;
  const PetscScalar *uv = NULL, *rv = NULL;
  PetscMPIInt      rank, size;
  PetscReal        local[4] = {0.0, 0.0, 0.0, 0.0}, global[4];
  const PetscReal  qv[2] = {-0.57735026918962576451, 0.57735026918962576451};
  PetscReal        jac = ctx->h * ctx->h * ctx->h / 8.0, dmap = 2.0 / ctx->h;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCall(VecScatterCreateToAll(u, &su, &ua));
  PetscCall(VecScatterCreateToAll(ref, &sr, &ra));
  PetscCall(VecScatterBegin(su, u, ua, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterEnd(su, u, ua, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterBegin(sr, ref, ra, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterEnd(sr, ref, ra, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecGetArrayRead(ua, &uv));
  PetscCall(VecGetArrayRead(ra, &rv));
  for (PetscInt e = rank; e < ctx->ne_x * ctx->ne_y * ctx->ne_z; e += size) {
    PetscInt ex = e % ctx->ne_x, tmp = e / ctx->ne_x;
    PetscInt ez = tmp % ctx->ne_z, ey = tmp / ctx->ne_z;
    PetscInt gids[8] = {
      NodeId(ctx, ex,     ey,     ez), NodeId(ctx, ex + 1, ey,     ez),
      NodeId(ctx, ex,     ey + 1, ez), NodeId(ctx, ex + 1, ey + 1, ez),
      NodeId(ctx, ex,     ey,     ez + 1), NodeId(ctx, ex + 1, ey,     ez + 1),
      NodeId(ctx, ex,     ey + 1, ez + 1), NodeId(ctx, ex + 1, ey + 1, ez + 1)
    };
    for (PetscInt a = 0; a < 2; a++) {
      for (PetscInt bb = 0; bb < 2; bb++) {
        for (PetscInt c = 0; c < 2; c++) {
          PetscReal N[8], dx[8], dy[8], dz[8];
          PetscScalar uq = 0.0, rq = 0.0, gu[3] = {0.0, 0.0, 0.0}, gr[3] = {0.0, 0.0, 0.0};
          Q1Shape(qv[a], qv[bb], qv[c], N, dx, dy, dz);
          for (PetscInt i = 0; i < 8; i++) {
            PetscScalar du = uv[gids[i]] - rv[gids[i]];
            uq += N[i] * du;
            rq += N[i] * rv[gids[i]];
            gu[0] += dmap * dx[i] * du; gu[1] += dmap * dy[i] * du; gu[2] += dmap * dz[i] * du;
            gr[0] += dmap * dx[i] * rv[gids[i]]; gr[1] += dmap * dy[i] * rv[gids[i]]; gr[2] += dmap * dz[i] * rv[gids[i]];
          }
          local[0] += jac * PetscAbsScalar(uq) * PetscAbsScalar(uq);
          local[1] += jac * PetscAbsScalar(rq) * PetscAbsScalar(rq);
          local[2] += jac * (ctx->omega * ctx->omega * PetscAbsScalar(uq) * PetscAbsScalar(uq) + PetscAbsScalar(gu[0]) * PetscAbsScalar(gu[0]) + PetscAbsScalar(gu[1]) * PetscAbsScalar(gu[1]) + PetscAbsScalar(gu[2]) * PetscAbsScalar(gu[2]));
          local[3] += jac * (ctx->omega * ctx->omega * PetscAbsScalar(rq) * PetscAbsScalar(rq) + PetscAbsScalar(gr[0]) * PetscAbsScalar(gr[0]) + PetscAbsScalar(gr[1]) * PetscAbsScalar(gr[1]) + PetscAbsScalar(gr[2]) * PetscAbsScalar(gr[2]));
        }
      }
    }
  }
  PetscCall(VecRestoreArrayRead(ra, &rv));
  PetscCall(VecRestoreArrayRead(ua, &uv));
  PetscCallMPI(MPI_Allreduce(local, global, 4, MPIU_REAL, MPI_SUM, PETSC_COMM_WORLD));
  *l2rel = PetscSqrtReal(global[0] / global[1]);
  *hk1rel = PetscSqrtReal(global[2] / global[3]);
  PetscCall(VecScatterDestroy(&sr));
  PetscCall(VecScatterDestroy(&su));
  PetscCall(VecDestroy(&ra));
  PetscCall(VecDestroy(&ua));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3Destroy(W3PC *w)
{
  PetscFunctionBeginUser;
  for (PetscInt j = 0; j < w->nlocal_patches; j++) PetscCall(W3PatchDestroy(&w->patches[j]));
  PetscCall(PetscFree(w->patches));
  PetscCall(VecDestroy(&w->coarse_corr));
  PetscCall(VecDestroy(&w->residual));
  PetscCall(VecDestroy(&w->ras));
  PetscCall(VecDestroy(&w->coarse_sol));
  PetscCall(VecDestroy(&w->coarse_rhs));
  PetscCall(W3CoarseMumpsSingleDestroy(&w->coarse_single));
  PetscCall(KSPDestroy(&w->coarse_ksp));
  PetscCall(MatDestroy(&w->Kc));
  PetscCall(MatDestroy(&w->APhi));
  PetscCall(MatDestroy(&w->Phi));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3SaveResult(Vec x, const char filename[])
{
  PetscFunctionBeginUser;
  if (filename[0]) PetscCall(SaveVecBinary(x, filename));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode W3FineDirect(const AppCtx *ctx, Mat A, Vec b, const char filename[])
{
  KSP            ksp = NULL;
  PC             pc;
  Vec            x = NULL, r = NULL, dx = NULL, xbest = NULL;
  PetscReal      bn, rn, best_rn = PETSC_MAX_REAL;
  PetscLogDouble t0, t1, t2, peak, peak_max = 0.0, peak_sum = 0.0;
  PetscInt       refine_steps = 0;

  PetscFunctionBeginUser;
  PetscCall(VecDuplicate(b, &x));
  PetscCall(VecDuplicate(b, &r));
  PetscCall(VecDuplicate(b, &dx));
  PetscCall(VecDuplicate(b, &xbest));
  PetscCall(KSPCreate(PETSC_COMM_WORLD, &ksp));
  PetscCall(KSPSetOperators(ksp, A, A));
  PetscCall(KSPSetType(ksp, KSPPREONLY));
  PetscCall(KSPGetPC(ksp, &pc));
  PetscCall(PCSetType(pc, PCLU));
  PetscCall(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
  PetscCall(KSPSetOptionsPrefix(ksp, "fine_direct_"));
  PetscCall(KSPSetFromOptions(ksp));
  PetscCall(PetscTime(&t0));
  PetscCall(KSPSetUp(ksp));
  PetscCall(PetscTime(&t1));
  PetscCall(KSPSolve(ksp, b, x));
  PetscCall(VecNorm(b, NORM_2, &bn));
  for (PetscInt it = 0; it <= 10; it++) {
    PetscCall(MatMult(A, x, r));
    PetscCall(VecAYPX(r, -1.0, b));
    PetscCall(VecNorm(r, NORM_2, &rn));
    if (rn < best_rn) {
      best_rn = rn;
      PetscCall(VecCopy(x, xbest));
    }
    if (rn / bn <= 1e-7) break;
    if (it == 10) break;
    PetscCall(KSPSolve(ksp, r, dx));
    PetscCall(VecAXPY(x, 1.0, dx));
    refine_steps++;
  }
  PetscCall(VecCopy(xbest, x));
  PetscCall(PetscTime(&t2));
  PetscCall(MatMult(A, x, r));
  PetscCall(VecAXPY(r, -1.0, b));
  PetscCall(VecNorm(r, NORM_2, &rn));
  PetscCall(W3SaveResult(x, filename));
  PetscCall(PetscMemoryGetMaximumUsage(&peak));
  PetscCallMPI(MPI_Reduce(&peak, &peak_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&peak, &peak_sum, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "fine direct MUMPS: |Ax-b|/|b| %.6e, refinement steps %" PetscInt_FMT ", factor(s) %.3f, solve+refine(s) %.3f, peak memory rank0/max-rank/sum-ranks %.3f/%.3f/%.3f GiB\n",
    (double)(rn / bn), refine_steps, (double)(t1 - t0), (double)(t2 - t1),
    (double)(peak / PetscPowReal(1024.0, 3.0)), (double)(peak_max / PetscPowReal(1024.0, 3.0)), (double)(peak_sum / PetscPowReal(1024.0, 3.0))));
  if (filename[0]) PetscCall(PetscPrintf(PETSC_COMM_WORLD, "saved fine reference: %s\n", filename));
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&xbest));
  PetscCall(VecDestroy(&dx));
  PetscCall(VecDestroy(&r));
  PetscCall(VecDestroy(&x));
  PetscFunctionReturn(PETSC_SUCCESS);
}

#include "continuous_outer.inc"

int main(int argc, char **argv)
{
  AppCtx         ctx;
  W3PC           w;
  Mat            A = NULL;
  Vec            b = NULL, x = NULL, xcoarse = NULL, r = NULL, ref = NULL;
  KSP            outer = NULL;
  PC             pc;
  PetscMPIInt    rank, size;
  PetscInt       fine_n = 64, coarse_n = 8, overlap = 2, dd_overlap = -1, basis_overlap = -1;
  PetscInt       basis_solve_overlap = -1, basis_support_overlap = -1, level = 0;
  PetscInt       outer_restart = 20, outer_max_it = 100;
  PetscReal      ppw = 20.0, outer_rtol = 1e-6, source_x = 0.5, source_y = 0.5, source_z = 0.5;
  PetscBool      fine_direct_only = PETSC_FALSE, coarse_only = PETSC_FALSE, fast_blended = PETSC_FALSE, blended_ibc_face = PETSC_FALSE;
  PetscBool      mixed_direct_precision = PETSC_TRUE;
  PetscBool      compare_outers = PETSC_FALSE;
  char           reference_file[PETSC_MAX_PATH_LEN] = "";
  char           solution_file[PETSC_MAX_PATH_LEN] = "";
  char           coarse_solution_file[PETSC_MAX_PATH_LEN] = "";
  PetscLogDouble t0, tA, tS, tC, tF, t1, t2, t3;
  PetscLogDouble local_setup, galerkin_setup, coarse_factor, coarse_solve = 0.0;
  PetscReal      coarse_relres = 0.0, bn, rn, coarse_l2 = -1.0, coarse_h1 = -1.0, outer_l2 = -1.0, outer_h1 = -1.0;
  PetscLogDouble peak, peak_max = 0.0, peak_sum = 0.0, local_pc_max = 0.0, coarse_pc_max = 0.0;
  PetscInt       its, ndof;
  KSPConvergedReason reason;

  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));
  PetscCall(PetscMemorySetGetMaximumUsage());
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));

  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "3D level-0 IBC-face wavelet EMs-HS", "");
  PetscCall(PetscOptionsInt("-fine_n", "fine elements in each direction", "", fine_n, &fine_n, NULL));
  PetscCall(PetscOptionsInt("-coarse_n", "coarse cells in each direction", "", coarse_n, &coarse_n, NULL));
  PetscCall(PetscOptionsInt("-overlap", "common overlap width used unless a separate DD or basis width is set", "", overlap, &overlap, NULL));
  PetscCall(PetscOptionsInt("-dd_overlap", "one-level RAS overlap width delta_DD/h in fine cells", "", dd_overlap, &dd_overlap, NULL));
  PetscCall(PetscOptionsInt("-basis_overlap", "backward-compatible common basis solve/support overlap", "", basis_overlap, &basis_overlap, NULL));
  PetscCall(PetscOptionsInt("-basis_solve_overlap", "basis local-solve oversampling width in fine cells", "", basis_solve_overlap, &basis_solve_overlap, NULL));
  PetscCall(PetscOptionsInt("-basis_support_overlap", "basis restriction and POU overlap width in fine cells", "", basis_support_overlap, &basis_support_overlap, NULL));
  PetscCall(PetscOptionsInt("-level", "wavelet level; this 3D implementation supports only ell=0", "", level, &level, NULL));
  PetscCall(PetscOptionsReal("-ppw", "fine-grid points per wavelength", "", ppw, &ppw, NULL));
  PetscCall(PetscOptionsReal("-source_x", "point-source x coordinate", "", source_x, &source_x, NULL));
  PetscCall(PetscOptionsReal("-source_y", "point-source y coordinate", "", source_y, &source_y, NULL));
  PetscCall(PetscOptionsReal("-source_z", "point-source z coordinate", "", source_z, &source_z, NULL));
  PetscCall(PetscOptionsInt("-outer_restart", "right-FGMRES restart", "", outer_restart, &outer_restart, NULL));
  PetscCall(PetscOptionsInt("-outer_max_it", "right-FGMRES maximum PC calls", "", outer_max_it, &outer_max_it, NULL));
  PetscCall(PetscOptionsReal("-outer_rtol", "true relative residual tolerance", "", outer_rtol, &outer_rtol, NULL));
  PetscCall(PetscOptionsBool("-fine_direct_only", "solve only the fine problem with MUMPS", "", fine_direct_only, &fine_direct_only, NULL));
  PetscCall(PetscOptionsBool("-coarse_only", "construct and solve only the coarse problem", "", coarse_only, &coarse_only, NULL));
  PetscCall(PetscOptionsBool("-fast_blended", "use the algebraically equivalent tensor-product blended element formula", "", fast_blended, &fast_blended, NULL));
  PetscCall(PetscOptionsBool("-blended_ibc_face", "use 0.5 GLL + 0.5 Gauss for the artificial IBC face mass", "", blended_ibc_face, &blended_ibc_face, NULL));
  PetscCall(PetscOptionsBool("-mixed_direct_precision", "use single-precision MUMPS for local and coarse direct solves", "", mixed_direct_precision, &mixed_direct_precision, NULL));
  PetscCall(PetscOptionsBool("-compare_outers", "reuse setup for physical BiCGStab and stationary iteration", "", compare_outers, &compare_outers, NULL));
  PetscCall(PetscOptionsString("-reference_file", "PETSc binary fine reference Vec", "", reference_file, reference_file, sizeof(reference_file), NULL));
  PetscCall(PetscOptionsString("-solution_file", "save EMs-HS solution Vec", "", solution_file, solution_file, sizeof(solution_file), NULL));
  PetscCall(PetscOptionsString("-coarse_solution_file", "save standalone coarse solution Vec", "", coarse_solution_file, coarse_solution_file, sizeof(coarse_solution_file), NULL));
  PetscOptionsEnd();

  if (dd_overlap < 0) dd_overlap = overlap;
  if (basis_overlap < 0) basis_overlap = overlap;
  if (basis_solve_overlap < 0) basis_solve_overlap = basis_overlap;
  if (basis_support_overlap < 0) basis_support_overlap = basis_overlap;

  PetscCheck(level == 0, PETSC_COMM_WORLD, PETSC_ERR_SUP, "3D study is intentionally restricted to wavelet level ell=0");
  PetscCheck(fine_n > 0 && coarse_n > 0 && fine_n % coarse_n == 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_INCOMP,
             "fine_n must be divisible by coarse_n");
  PetscCheck(dd_overlap > 0 && dd_overlap <= fine_n / coarse_n, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
             "dd_overlap must satisfy 1 <= delta_DD/h <= H/h");
  PetscCheck(basis_solve_overlap > 0 && basis_solve_overlap <= fine_n / coarse_n, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
             "basis_solve_overlap must satisfy 1 <= delta_solve/h <= H/h");
  PetscCheck(basis_support_overlap > 0 && basis_support_overlap <= basis_solve_overlap, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
             "basis_support_overlap must satisfy 1 <= delta_support/h <= delta_solve/h");
  PetscCheck(ppw > 0.0 && outer_rtol > 0.0 && outer_restart > 0 && outer_max_it > 0, PETSC_COMM_WORLD, PETSC_ERR_ARG_OUTOFRANGE,
             "bad ppw or outer solver parameters");

  PetscCall(PetscMemzero(&ctx, sizeof(ctx)));
  ctx.nphys = ctx.nex = ctx.ney = ctx.nez = fine_n;
  ctx.ne_x = ctx.ne_y = ctx.ne_z = fine_n;
  ctx.nx_opt = ctx.ny_opt = ctx.nz_opt = fine_n;
  ctx.npml = ctx.npmly = 0;
  ctx.h = 1.0 / (PetscReal)fine_n;
  ctx.h_m = 0.0;
  ctx.v0 = ctx.vel_min = ctx.vel_max = 1.0;
  ctx.vel_h = ctx.h;
  ctx.min_ppw = ppw;
  ctx.frequency = 1.0 / (ctx.h * ppw);
  ctx.omega = 2.0 * PETSC_PI * ctx.frequency;
  ctx.operator_scale = 1.0;
  ctx.mass_shift_beta = 0.0;
  ctx.abc_sign = 1.0;
  ctx.use_1st_abc = PETSC_TRUE;
  ctx.fast_blended_element = fast_blended;
  ctx.gauss_quadrature = PETSC_TRUE;
  ctx.source_point_delta = PETSC_TRUE;
  ctx.source_set = PETSC_TRUE;
  ctx.source_x_m = source_x;
  ctx.source_y_m = source_y;
  ctx.source_z_m = source_z;
  ctx.source_support = ctx.h;
  ndof = (fine_n + 1) * (fine_n + 1) * (fine_n + 1);

  if (fine_direct_only) {
    PetscCall(HmsSetOptionIfMissing("-fine_direct_mat_mumps_icntl_14", "160"));
    PetscCall(HmsSetOptionIfMissing("-fine_direct_mat_mumps_icntl_29", "2"));
  } else {
    if (mixed_direct_precision) {
      PetscCall(HmsSetOptionIfMissing("-wavelet_local_pc_precision", "single"));
    }
    PetscCall(HmsSetOptionIfMissing("-wavelet_local_mat_mumps_icntl_14", "160"));
    if (!mixed_direct_precision) {
      PetscCall(HmsSetOptionIfMissing("-wavelet_coarse_mat_mumps_icntl_14", "160"));
      PetscCall(HmsSetOptionIfMissing("-wavelet_coarse_mat_mumps_icntl_29", "2"));
    }
  }

  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "3D Helmholtz test: fine=%" PetscInt_FMT "^3, coarse=%" PetscInt_FMT "^3, H/h=%" PetscInt_FMT ", ell=0, DD delta/h=%" PetscInt_FMT ", basis solve/support delta/h=%" PetscInt_FMT "/%" PetscInt_FMT ", fine ppw=%.6g, coarse ppw=%.6g, k=%.6g, MPI=%d, quadrature=Gauss(2)^3/Gauss(2)^2, artificial IBC face=%s, physical/local signs=+ik/+iH^{-1}\n",
    fine_n, coarse_n, fine_n / coarse_n, dd_overlap, basis_solve_overlap, basis_support_overlap, (double)ppw,
    (double)(ppw / ((PetscReal)fine_n / (PetscReal)coarse_n)), (double)ctx.omega, (int)size,
    blended_ibc_face ? "blended" : "Gauss"));
  if (compare_outers) PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "outer comparison: FGMRES(20) and stationary, direct |Ax-b|/|b| <= %.6g\n", (double)outer_rtol));
  else PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "outer solver: right-preconditioned FGMRES(%" PetscInt_FMT "), stopping on direct |Ax-b|/|b| <= %.6g\n",
    outer_restart, (double)outer_rtol));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "direct-solve precision: local/coarse MUMPS=%s, outer operator and Krylov basis=double\n",
    mixed_direct_precision ? "single" : "double"));

  PetscCall(PetscTime(&t0));
  PetscCall(AssembleSystem(&ctx, &A, &b));
  PetscCall(PetscTime(&tA));
  if (fine_direct_only) {
    PetscCall(W3FineDirect(&ctx, A, b, reference_file));
    PetscCall(VecDestroy(&b));
    PetscCall(MatDestroy(&A));
    PetscCall(PetscFinalize());
    return 0;
  }

  if (reference_file[0]) PetscCall(LoadVecBinary(PETSC_COMM_WORLD, ndof, reference_file, &ref));
  PetscCall(W3Setup(&ctx, A, b, coarse_n, dd_overlap, basis_solve_overlap, basis_support_overlap,
                    blended_ibc_face, mixed_direct_precision,
                    &w, &local_setup, &galerkin_setup, &coarse_factor));
  PetscCall(PetscTime(&tS));
  PetscCall(VecDuplicate(b, &xcoarse));
  PetscCall(W3CoarseSolution(&w, b, xcoarse, &coarse_relres, &coarse_solve));
  if (ref) PetscCall(W3RelativeFeErrors(&ctx, xcoarse, ref, &coarse_l2, &coarse_h1));
  PetscCall(W3SaveResult(xcoarse, coarse_solution_file));
  PetscCall(PetscTime(&tC));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "standalone coarse solution: |Ax0-b|/|b| %.6e, relative L2 %.6e, relative Hk1 %.6e, solve(s) %.3f\n",
    (double)coarse_relres, (double)coarse_l2, (double)coarse_h1, (double)coarse_solve));

  if (!coarse_only && compare_outers) {
    Vec dx;
    PetscCall(VecDuplicate(b, &x));
    PetscCall(VecDuplicate(b, &r));
    PetscCall(VecDuplicate(b, &dx));
    PetscCall(VecNorm(b, NORM_2, &bn));
    w.pc_calls = 0;
    PetscCall(PetscTime(&t1));
    PetscCall(W3FGMRES20(A, b, x, &w, outer_max_it, outer_rtol, &its, &rn, &reason));
    PetscCall(PetscTime(&t2));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "OUTER_RESULT fgmres20 N=%" PetscInt_FMT " PC=%" PetscInt_FMT " reason=%d residual=%.12e solve=%.3f\n",
      fine_n, w.pc_calls, (int)reason, (double)(rn / PetscMax(bn, PETSC_SMALL)), (double)(t2-t1)));
    PetscCall(PCCreate(PETSC_COMM_WORLD, &pc));
    PetscCall(PCSetOperators(pc, A, A));
    PetscCall(PCSetType(pc, PCSHELL));
    PetscCall(PCShellSetContext(pc, &w));
    PetscCall(PCShellSetApply(pc, W3PCApply));
    PetscCall(PCSetUp(pc));
    PetscCall(VecSet(x, 0.0));
    w.pc_calls = 0;
    reason = KSP_DIVERGED_ITS;
    PetscCall(PetscTime(&t1));
    for (;;) {
      PetscCall(MatMult(A, x, r));
      PetscCall(VecAYPX(r, -1.0, b));
      PetscCall(VecNorm(r, NORM_2, &rn));
      PetscCall(PetscPrintf(PETSC_COMM_WORLD, "STATIONARY PC=%" PetscInt_FMT " residual=%.12e\n",
        w.pc_calls, (double)(rn / PetscMax(bn, PETSC_SMALL))));
      if (rn <= outer_rtol*bn) { reason=KSP_CONVERGED_RTOL; break; }
      if (PetscIsInfOrNanReal(rn) || rn > 1e8*bn) { reason=KSP_DIVERGED_DTOL; break; }
      PetscCall(PetscTime(&t2));
      if (w.pc_calls >= outer_max_it || t2-t1 > 600.0) break;
      PetscCall(PCApply(pc, r, dx));
      PetscCall(VecAXPY(x, 1.0, dx));
    }
    PetscCall(PetscTime(&t3));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "OUTER_RESULT stationary N=%" PetscInt_FMT " PC=%" PetscInt_FMT " reason=%d residual=%.12e solve=%.3f\n",
      fine_n, w.pc_calls, (int)reason, (double)(rn / PetscMax(bn, PETSC_SMALL)), (double)(t3-t1)));
    PetscCall(PCDestroy(&pc));
    PetscCall(VecDestroy(&dx));
  } else if (!coarse_only) {
    W3TrueResidualCtx true_conv = {0};
    PetscInt          residual_restarts = 0, total_its = 0;
    PetscBool         outer_converged = PETSC_FALSE;

    PetscCall(VecDuplicate(b, &x));
    PetscCall(VecDuplicate(b, &r));
    PetscCall(VecDuplicate(b, &true_conv.solution));
    PetscCall(VecDuplicate(b, &true_conv.residual));
    PetscCall(VecSet(true_conv.solution, 0.0));
    PetscCall(VecNorm(b, NORM_2, &bn));
    true_conv.target_abs = outer_rtol * bn;
    true_conv.last_abs = PETSC_MAX_REAL;
    true_conv.max_it = outer_max_it;
    PetscCall(KSPCreate(PETSC_COMM_WORLD, &outer));
    PetscCall(KSPSetOperators(outer, A, A));
    PetscCall(KSPSetType(outer, KSPFGMRES));
    PetscCall(KSPGMRESSetRestart(outer, outer_restart));
    PetscCall(KSPSetTolerances(outer, outer_rtol, PETSC_DEFAULT, PETSC_DEFAULT, outer_max_it));
    PetscCall(KSPSetPCSide(outer, PC_RIGHT));
    PetscCall(KSPSetNormType(outer, KSP_NORM_UNPRECONDITIONED));
    PetscCall(KSPSetOptionsPrefix(outer, "outer_"));
    PetscCall(KSPGetPC(outer, &pc));
    PetscCall(PCSetType(pc, PCSHELL));
    PetscCall(PCShellSetContext(pc, &w));
    PetscCall(PCShellSetApply(pc, W3PCApply));
    PetscCall(KSPSetFromOptions(outer));
    PetscCall(KSPSetConvergenceTest(outer, W3TrueResidualConverged, &true_conv, NULL));
    PetscCall(VecSet(x, 0.0));
    PetscCall(PetscTime(&t1));
    for (;;) {
      PetscInt remaining = outer_max_it - w.pc_calls;

      if (remaining <= 0) {
        reason = KSP_DIVERGED_ITS;
        break;
      }
      PetscCall(KSPSetInitialGuessNonzero(outer, residual_restarts ? PETSC_TRUE : PETSC_FALSE));
      PetscCall(KSPSetTolerances(outer, outer_rtol, PETSC_DEFAULT, PETSC_DEFAULT, remaining));
      true_conv.last_abs = PETSC_MAX_REAL;
      true_conv.max_it = remaining;
      PetscCall(KSPSolve(outer, b, x));
      PetscCall(KSPGetIterationNumber(outer, &its));
      PetscCall(KSPGetConvergedReason(outer, &reason));
      total_its += its;
      PetscCall(MatMult(A, x, r));
      PetscCall(VecAXPY(r, -1.0, b));
      PetscCall(VecNorm(r, NORM_2, &rn));
      if (rn <= outer_rtol * bn) {
        outer_converged = PETSC_TRUE;
        break;
      }
      if (reason < 0 || residual_restarts >= 8) break;
      residual_restarts++;
    }
    PetscCall(PetscTime(&t2));
    its = total_its;
    if (ref) PetscCall(W3RelativeFeErrors(&ctx, x, ref, &outer_l2, &outer_h1));
    PetscCall(W3SaveResult(x, solution_file));
    PetscCall(PetscTime(&t3));
    PetscCallMPI(MPI_Reduce(&w.pc_local_solve_time, &local_pc_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Reduce(&w.pc_coarse_solve_time, &coarse_pc_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "EMs-HS: converged %s, PC calls %" PetscInt_FMT ", FGMRES its %" PetscInt_FMT ", residual restarts %" PetscInt_FMT ", reason %d, |Ax-b|/|b| %.6e, relative L2 %.6e, relative Hk1 %.6e\n",
      outer_converged ? "yes" : "no", w.pc_calls, its, residual_restarts, (int)reason,
      (double)(rn / bn), (double)outer_l2, (double)outer_h1));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
      "solve(s) %.3f, check/save(s) %.3f, PC local max-rank total(s) %.3f, PC coarse/residual max-rank total(s) %.3f\n",
      (double)(t2 - t1), (double)(t3 - t2), (double)local_pc_max, (double)coarse_pc_max));
    PetscCall(VecDestroy(&true_conv.residual));
    PetscCall(VecDestroy(&true_conv.solution));
  } else {
    t3 = tC;
    its = 0;
  }
  PetscCall(PetscTime(&tF));
  PetscCall(PetscMemoryGetMaximumUsage(&peak));
  PetscCallMPI(MPI_Reduce(&peak, &peak_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&peak, &peak_sum, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "setup(s): fine %.3f, local basis+MUMPS %.3f, Galerkin %.3f, coarse MUMPS %.3f, total %.3f; coarse postcheck %.3f; final bookkeeping %.3f\n",
    (double)(tA - t0), (double)local_setup, (double)galerkin_setup, (double)coarse_factor,
    (double)(tS - t0), (double)(tC - tS), (double)(tF - t3)));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
    "peak memory rank0/max-rank/sum-ranks %.3f/%.3f/%.3f GiB\n",
    (double)(peak / PetscPowReal(1024.0, 3.0)), (double)(peak_max / PetscPowReal(1024.0, 3.0)), (double)(peak_sum / PetscPowReal(1024.0, 3.0))));

  PetscCall(KSPDestroy(&outer));
  PetscCall(VecDestroy(&ref));
  PetscCall(VecDestroy(&r));
  PetscCall(VecDestroy(&x));
  PetscCall(VecDestroy(&xcoarse));
  PetscCall(W3Destroy(&w));
  PetscCall(VecDestroy(&b));
  PetscCall(MatDestroy(&A));
  PetscCall(PetscFinalize());
  return 0;
}
