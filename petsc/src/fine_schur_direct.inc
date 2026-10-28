/* Exact one-level substructuring for the fine Q1 problem.
   Coarse cells define nonoverlapping substructures. Their fine-interior
   unknowns are eliminated by local PARDISO solves; the assembled interface
   Schur complement is solved by the distributed MUMPS path. */

static PetscErrorCode FineSchurBuildLocalMaps(const MsCtx *ms, PetscInt **imap_out, PetscInt **bmap_out, PetscInt **blids_out)
{
  PetscInt *imap = NULL, *bmap = NULL, *blids = NULL;
  PetscInt  ni = 0, nb = 0;

  PetscFunctionBeginUser;
  PetscCall(PetscMalloc1(ms->local_fine_nodes, &imap));
  PetscCall(PetscMalloc1(ms->local_fine_nodes, &bmap));
  PetscCall(PetscMalloc1(ms->local_p_bdry, &blids));
  for (PetscInt iy = 0; iy <= ms->by; iy++) {
    for (PetscInt iz = 0; iz <= ms->bz; iz++) {
      for (PetscInt ix = 0; ix <= ms->bx; ix++) {
        PetscInt lid = LocalFineNodeId(ms, ix, iy, iz);
        if (IsBoundaryFineNode(ms, ix, iy, iz)) {
          imap[lid] = -1;
          bmap[lid] = nb;
          blids[nb++] = lid;
        } else {
          imap[lid] = ni++;
          bmap[lid] = -1;
        }
      }
    }
  }
  PetscCheck(ni == ms->local_p_inner && nb == ms->local_p_bdry, PETSC_COMM_SELF, PETSC_ERR_PLIB,
             "Fine Schur local partition mismatch: interior %" PetscInt_FMT "/%" PetscInt_FMT ", boundary %" PetscInt_FMT "/%" PetscInt_FMT,
             ni, ms->local_p_inner, nb, ms->local_p_bdry);
  *imap_out  = imap;
  *bmap_out  = bmap;
  *blids_out = blids;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FineSchurPardisoSolve(const MsCtx *ms, PetscInt n_, PetscInt nrhs_, int ia[], int ja[], PetscScalar aval[], PetscScalar rhs[], PetscScalar sol[], PetscLogDouble *factor_solve_time)
{
  PardisoFunc pardiso = GetPardisoFunc();
  void       *pt[64];
  int         iparm[64], maxfct = 1, mnum = 1, mtype = 13, phase = 13;
  int         n = (int)n_, nrhs = (int)nrhs_, msglvl = 0, error = 0;
  int        *perm = NULL;
  PetscLogDouble t0, t1;

  PetscFunctionBeginUser;
  PetscCheck(pardiso, PETSC_COMM_SELF, PETSC_ERR_SUP, "Fine Schur local PARDISO symbol is unavailable: %s",
             hms_pardiso_lookup_error[0] ? hms_pardiso_lookup_error : "none");
  PetscCheck(n_ > 0 && n_ <= INT_MAX && nrhs_ > 0 && nrhs_ <= INT_MAX, PETSC_COMM_SELF, PETSC_ERR_SUP,
             "Fine Schur local PARDISO dimensions exceed 32-bit local indices");
  ConfigureLocalInteriorPardiso(ms, iparm, pt);
  PetscCall(PetscCalloc1(n_, &perm));
  PetscCall(PetscTime(&t0));
  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, aval, ia, ja, perm, &nrhs, iparm, &msglvl, rhs, sol, &error);
  PetscCall(PetscTime(&t1));
  PetscCheck(error == 0, PETSC_COMM_SELF, PETSC_ERR_LIB,
             "Fine Schur local PARDISO factor/solve failed with error=%d, n=%d, nrhs=%d", error, n, nrhs);
  phase = -1;
  pardiso(pt, &maxfct, &mnum, &mtype, &phase, &n, aval, ia, ja, perm, &nrhs, iparm, &msglvl, rhs, sol, &error);
  PetscCheck(error == 0, PETSC_COMM_SELF, PETSC_ERR_LIB, "Fine Schur local PARDISO release failed with error=%d", error);
  if (factor_solve_time) *factor_solve_time += t1 - t0;
  PetscCall(PetscFree(perm));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FineSchurBuildLocalBlocks(const AppCtx *ctx, const MsCtx *ms, PetscInt cx, PetscInt cy, PetscInt cz,
                                                PetscInt **imap_out, PetscInt **bmap_out, PetscInt **blids_out,
                                                int **ia_out, int **ja_out, PetscScalar **aii_out,
                                                PetscScalar **aib_col_out, PetscScalar **abi_col_out,
                                                PetscScalar **abb_col_out, PetscScalar **f0_out)
{
  const PetscInt cap = 27;
  PetscInt       ni = ms->local_p_inner, nb = ms->local_p_bdry;
  PetscInt      *imap = NULL, *bmap = NULL, *blids = NULL;
  int           *counts = NULL, *cols = NULL, *ia = NULL, *ja = NULL;
  PetscScalar   *vals = NULL, *aii = NULL, *aib = NULL, *abi = NULL, *abb = NULL, *f0 = NULL;
  PetscInt       nnz = 0, pos = 0;

  PetscFunctionBeginUser;
  PetscCall(FineSchurBuildLocalMaps(ms, &imap, &bmap, &blids));
  PetscCall(PetscCalloc1(ni, &counts));
  PetscCall(PetscMalloc1((size_t)ni * (size_t)cap, &cols));
  PetscCall(PetscCalloc1((size_t)ni * (size_t)cap, &vals));
  PetscCall(PetscCalloc1((size_t)ni * (size_t)nb, &aib));
  PetscCall(PetscCalloc1((size_t)nb * (size_t)ni, &abi));
  PetscCall(PetscCalloc1((size_t)nb * (size_t)nb, &abb));

  for (PetscInt ly = 0; ly < ms->by; ly++) {
    for (PetscInt lz = 0; lz < ms->bz; lz++) {
      for (PetscInt lx = 0; lx < ms->bx; lx++) {
        PetscInt    ex = cx * ms->bx + lx, ey = cy * ms->by + ly, ez = cz * ms->bz + lz;
        PetscInt    lids[8];
        PetscScalar Ke[64];

        lids[0] = LocalFineNodeId(ms, lx,     ly,     lz);
        lids[1] = LocalFineNodeId(ms, lx + 1, ly,     lz);
        lids[2] = LocalFineNodeId(ms, lx,     ly + 1, lz);
        lids[3] = LocalFineNodeId(ms, lx + 1, ly + 1, lz);
        lids[4] = LocalFineNodeId(ms, lx,     ly,     lz + 1);
        lids[5] = LocalFineNodeId(ms, lx + 1, ly,     lz + 1);
        lids[6] = LocalFineNodeId(ms, lx,     ly + 1, lz + 1);
        lids[7] = LocalFineNodeId(ms, lx + 1, ly + 1, lz + 1);
        PetscCall(ComputeElementMatrixOnly(ctx, ex, ey, ez, Ke));
        for (PetscInt a = 0; a < 8; a++) {
          PetscInt ri = imap[lids[a]], bi = bmap[lids[a]];
          for (PetscInt b = 0; b < 8; b++) {
            PetscInt    rj = imap[lids[b]], bj = bmap[lids[b]];
            PetscScalar kab = Ke[8 * a + b];
            if (ri >= 0) {
              if (rj >= 0) PetscCall(AddLocalCsrEntry(ri, rj, kab, cap, counts, cols, vals));
              else aib[(size_t)bj * (size_t)ni + (size_t)ri] += kab;
            } else {
              if (rj >= 0) abi[(size_t)rj * (size_t)nb + (size_t)bi] += kab;
              else abb[(size_t)bj * (size_t)nb + (size_t)bi] += kab;
            }
          }
        }
      }
    }
  }

  for (PetscInt row = 0; row < ni; row++) {
    for (int i = 1; i < counts[row]; i++) {
      int         ctmp = cols[(size_t)row * cap + i], j = i - 1;
      PetscScalar vtmp = vals[(size_t)row * cap + i];
      while (j >= 0 && cols[(size_t)row * cap + j] > ctmp) {
        cols[(size_t)row * cap + j + 1] = cols[(size_t)row * cap + j];
        vals[(size_t)row * cap + j + 1] = vals[(size_t)row * cap + j];
        j--;
      }
      cols[(size_t)row * cap + j + 1] = ctmp;
      vals[(size_t)row * cap + j + 1] = vtmp;
    }
    nnz += counts[row];
  }
  PetscCall(PetscMalloc1(ni + 1, &ia));
  PetscCall(PetscMalloc1(nnz, &ja));
  PetscCall(PetscMalloc1(nnz, &aii));
  ia[0] = 1;
  for (PetscInt row = 0; row < ni; row++) {
    for (int q = 0; q < counts[row]; q++) {
      ja[pos]  = cols[(size_t)row * cap + q] + 1;
      aii[pos] = vals[(size_t)row * cap + q];
      pos++;
    }
    ia[row + 1] = (int)pos + 1;
  }
  PetscCall(BuildLocalFineRHS(ctx, ms, cx, cy, cz, &f0));

  *imap_out = imap; *bmap_out = bmap; *blids_out = blids;
  *ia_out = ia; *ja_out = ja; *aii_out = aii;
  *aib_col_out = aib; *abi_col_out = abi; *abb_col_out = abb; *f0_out = f0;
  PetscCall(PetscFree(vals));
  PetscCall(PetscFree(cols));
  PetscCall(PetscFree(counts));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FineSchurEliminateCell(const AppCtx *ctx, const MsCtx *ms, PetscInt cx, PetscInt cy, PetscInt cz,
                                             PetscScalar **srow_out, PetscScalar **fshur_out, PetscLogDouble *local_solve_time)
{
  PetscInt     *imap = NULL, *bmap = NULL, *blids = NULL, *ilids = NULL;
  int          *ia = NULL, *ja = NULL;
  PetscScalar  *aii = NULL, *aib = NULL, *abi = NULL, *scol = NULL, *f0 = NULL;
  PetscScalar  *rhs = NULL, *sol = NULL, *srow = NULL, *fshur = NULL;
  PetscInt      ni = ms->local_p_inner, nb = ms->local_p_bdry, nrhs = nb + 1;

  PetscFunctionBeginUser;
  PetscCall(FineSchurBuildLocalBlocks(ctx, ms, cx, cy, cz, &imap, &bmap, &blids, &ia, &ja, &aii, &aib, &abi, &scol, &f0));
  PetscCall(PetscMalloc1(ni, &ilids));
  for (PetscInt lid = 0; lid < ms->local_fine_nodes; lid++) if (imap[lid] >= 0) ilids[imap[lid]] = lid;
  PetscCall(PetscMalloc1((size_t)ni * (size_t)nrhs, &rhs));
  PetscCall(PetscMalloc1((size_t)ni * (size_t)nrhs, &sol));
  PetscCall(PetscArraycpy(rhs, aib, (size_t)ni * (size_t)nb));
  for (PetscInt i = 0; i < ni; i++) rhs[(size_t)nb * (size_t)ni + (size_t)i] = f0[ilids[i]];
  PetscCall(FineSchurPardisoSolve(ms, ni, nrhs, ia, ja, aii, rhs, sol, local_solve_time));
  PetscCall(Hms2gGemm('N', 'N', nb, nb, ni, -1.0, abi, nb, sol, ni, 1.0, scol, nb));
  PetscCall(PetscMalloc1((size_t)nb * (size_t)nb, &srow));
  PetscCall(PetscMalloc1(nb, &fshur));
  for (PetscInt i = 0; i < nb; i++) {
    fshur[i] = f0[blids[i]];
    for (PetscInt k = 0; k < ni; k++) fshur[i] -= abi[(size_t)k * (size_t)nb + (size_t)i] * sol[(size_t)nb * (size_t)ni + (size_t)k];
    for (PetscInt j = 0; j < nb; j++) srow[(size_t)i * (size_t)nb + (size_t)j] = scol[(size_t)j * (size_t)nb + (size_t)i];
  }
  *srow_out = srow;
  *fshur_out = fshur;
  PetscCall(PetscFree(ilids));
  PetscCall(PetscFree(sol)); PetscCall(PetscFree(rhs));
  PetscCall(PetscFree(f0)); PetscCall(PetscFree(scol)); PetscCall(PetscFree(abi)); PetscCall(PetscFree(aib)); PetscCall(PetscFree(aii));
  PetscCall(PetscFree(ja)); PetscCall(PetscFree(ia));
  PetscCall(PetscFree(blids)); PetscCall(PetscFree(bmap)); PetscCall(PetscFree(imap));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode FineSchurRecoverCell(const AppCtx *ctx, const MsCtx *ms, PetscInt cx, PetscInt cy, PetscInt cz,
                                           const PetscScalar cb[], PetscScalar **ui_out, PetscLogDouble *local_solve_time)
{
  PetscInt     *imap = NULL, *bmap = NULL, *blids = NULL;
  int          *ia = NULL, *ja = NULL;
  PetscScalar  *aii = NULL, *aib = NULL, *abi = NULL, *abb = NULL, *f0 = NULL, *rhs = NULL, *ui = NULL;
  PetscInt      ni = ms->local_p_inner, nb = ms->local_p_bdry;

  PetscFunctionBeginUser;
  PetscCall(FineSchurBuildLocalBlocks(ctx, ms, cx, cy, cz, &imap, &bmap, &blids, &ia, &ja, &aii, &aib, &abi, &abb, &f0));
  PetscCall(PetscMalloc1(ni, &rhs));
  PetscCall(PetscMalloc1(ni, &ui));
  for (PetscInt lid = 0; lid < ms->local_fine_nodes; lid++) if (imap[lid] >= 0) rhs[imap[lid]] = f0[lid];
  for (PetscInt j = 0; j < nb; j++) {
    PetscScalar v = cb[j];
    for (PetscInt i = 0; i < ni; i++) rhs[i] -= aib[(size_t)j * (size_t)ni + (size_t)i] * v;
  }
  PetscCall(FineSchurPardisoSolve(ms, ni, 1, ia, ja, aii, rhs, ui, local_solve_time));
  *ui_out = ui;
  PetscCall(PetscFree(rhs)); PetscCall(PetscFree(f0)); PetscCall(PetscFree(abb)); PetscCall(PetscFree(abi)); PetscCall(PetscFree(aib)); PetscCall(PetscFree(aii));
  PetscCall(PetscFree(ja)); PetscCall(PetscFree(ia));
  PetscCall(PetscFree(blids)); PetscCall(PetscFree(bmap)); PetscCall(PetscFree(imap));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode SolveFineSchurDirect(const AppCtx *ctx, MsCtx *ms)
{
  PetscMPIInt    rank, size;
  PetscInt      *shur_map = NULL, *bdry = NULL, *inner = NULL, *gids = NULL;
  PetscInt       ncells, local_cells, cell_begin, cell_end, prealloc_nz;
  Mat            S = NULL, Afull = NULL;
  Vec            F = NULL, U = NULL, R = NULL, Uall = NULL, ufine = NULL, bfull = NULL, rfull = NULL;
  VecScatter     scat = NULL;
  KSP            coarse = NULL;
  PetscLogDouble tp0, tp1, t0, t1, t2, t3, t4, t5, t6, local_solve = 0.0, local_solve_max = 0.0, recover_solve = 0.0, recover_solve_max = 0.0;
  PetscReal      fnorm = 0.0, rnorm = 0.0, full_bnorm = 0.0, full_rnorm = -1.0, reference_rel = -1.0;
  PetscLogDouble peak = 0.0, peak_max = 0.0, peak_sum = 0.0;
  PetscBool      check_full = PETSC_FALSE;
  PetscBool      exact_prealloc = PETSC_TRUE;
  PetscBool      analysis_only = PETSC_FALSE;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  PetscCheck(ms->bx == ms->by && ms->bx == ms->bz, PETSC_COMM_WORLD, PETSC_ERR_SUP,
             "Fine exact-Schur prototype currently needs the same H/h in x,y,z");
  PetscCheck(!ms->top_dirichlet, PETSC_COMM_WORLD, PETSC_ERR_SUP,
             "Fine exact-Schur prototype currently targets full PML/ABC, not row-zeroed top Dirichlet");
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-fine_schur_check_full", &check_full, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-fine_schur_exact_prealloc", &exact_prealloc, NULL));
  PetscCall(PetscOptionsGetBool(NULL, NULL, "-fine_schur_analysis_only", &analysis_only, NULL));
  ms->pms = ms->bx;
  ms->local_p_nodes  = ms->local_fine_nodes;
  ms->local_p_inner  = (ms->bx - 1) * (ms->by - 1) * (ms->bz - 1);
  ms->local_p_bdry   = ms->local_p_nodes - ms->local_p_inner;
  ms->global_p_nodes = (ctx->ne_x + 1) * (ctx->ne_y + 1) * (ctx->ne_z + 1);
  ncells = ms->coarse_x * ms->coarse_y * ms->coarse_z;
  ms->shur_nodes = ms->global_p_nodes - ncells * ms->local_p_inner;
  prealloc_nz = PetscMin(ms->shur_nodes, 2 * ms->local_p_bdry);
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-fine_schur_prealloc_nz", &prealloc_nz, NULL));
  prealloc_nz = PetscMax((PetscInt)1, PetscMin(prealloc_nz, ms->shur_nodes));

  PetscCall(BuildGlobalShurMap(ms, &shur_map));
  PetscCall(PetscMalloc2(ms->local_p_bdry, &bdry, PetscMax(ms->local_p_inner, (PetscInt)1), &inner));
  PetscCall(PetscMalloc1(ms->local_p_bdry, &gids));
  PetscCall(BuildLocalPPartition(ms, bdry, inner));
  PetscCall(PetscTime(&tp0));
  if (exact_prealloc) {
    PetscCall(CreateCoarseAIJWithExactPattern(ms, shur_map, rank, size, &S));
  } else {
    PetscCall(MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, ms->shur_nodes, ms->shur_nodes,
                           prealloc_nz, NULL, prealloc_nz, NULL, &S));
    PetscCall(MatSetOption(S, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
  }
  PetscCall(PetscTime(&tp1));
  PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, ms->shur_nodes, &F));
  PetscCall(VecSet(F, 0.0));
  MsOwnedCellRange(ncells, rank, size, ms->ms_spatial_cell_partition, &cell_begin, &cell_end);
  local_cells = MsOwnedCellCount(ncells, rank, size, ms->ms_spatial_cell_partition);

  PetscCall(PetscTime(&t0));
  for (PetscInt li = 0; li < local_cells; li++) {
    PetscInt     c = MsOwnedCellId(li, rank, size, ms->ms_spatial_cell_partition, cell_begin);
    PetscInt     cx = c % ms->coarse_x, tmp = c / ms->coarse_x, cz = tmp % ms->coarse_z, cy = tmp / ms->coarse_z;
    PetscScalar *Sloc = NULL, *Floc = NULL;
    PetscCall(BuildCellBoundaryGids(ms, shur_map, bdry, cx, cy, cz, gids));
    PetscCall(FineSchurEliminateCell(ctx, ms, cx, cy, cz, &Sloc, &Floc, &local_solve));
    PetscCall(MatSetValues(S, ms->local_p_bdry, gids, ms->local_p_bdry, gids, Sloc, ADD_VALUES));
    PetscCall(VecSetValues(F, ms->local_p_bdry, gids, Floc, ADD_VALUES));
    PetscCall(PetscFree(Floc));
    PetscCall(PetscFree(Sloc));
  }
  PetscCall(PetscTime(&t1));
  PetscCall(MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY));
  PetscCall(VecAssemblyBegin(F));
  PetscCall(MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY));
  PetscCall(VecAssemblyEnd(F));
  PetscCall(PetscTime(&t2));
  PetscCallMPI(MPI_Reduce(&local_solve, &local_solve_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));

  if (analysis_only) {
    Mat           factor = NULL;
    MatFactorInfo info;
    PetscInt      est_max_mb = -1, est_sum_mb = -1, est_factor_entries = -1;

    PetscCall(HmsMaybeSetCoarseMumpsBlrOptions(ctx));
    PetscCall(MatGetFactor(S, MATSOLVERMUMPS, MAT_FACTOR_LU, &factor));
    PetscCall(MatSetOptionsPrefix(factor, "hms_coarse_"));
    PetscCall(MatFactorInfoInitialize(&info));
    PetscCall(PetscTime(&t3));
    PetscCall(MatLUFactorSymbolic(factor, S, NULL, NULL, &info));
    PetscCall(PetscTime(&t4));
    PetscCall(MatMumpsGetInfog(factor, 16, &est_max_mb));
    PetscCall(MatMumpsGetInfog(factor, 17, &est_sum_mb));
    PetscCall(MatMumpsGetInfog(factor, 20, &est_factor_entries));
    PetscCall(PetscMemoryGetMaximumUsage(&peak));
    PetscCallMPI(MPI_Reduce(&peak, &peak_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
    PetscCallMPI(MPI_Reduce(&peak, &peak_sum, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Fine exact-interface Schur MUMPS analysis only\n"));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
               "fine dofs %" PetscInt_FMT ", interface dofs %" PetscInt_FMT " (%.3f%%), coarse cells %" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT ", H/h %" PetscInt_FMT ", local interior/boundary %" PetscInt_FMT "/%" PetscInt_FMT "\n",
               ms->global_p_nodes, ms->shur_nodes, 100.0 * (double)ms->shur_nodes / (double)ms->global_p_nodes,
               ms->coarse_y, ms->coarse_z, ms->coarse_x, ms->bx, ms->local_p_inner, ms->local_p_bdry));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
               "time(s) interface pattern %.3f, local eliminate+Schur %.3f (PARDISO max-rank %.3f), interface MatAssembly %.3f, MUMPS analysis %.3f\n",
               (double)(tp1 - tp0), (double)(t1 - t0), (double)local_solve_max, (double)(t2 - t1), (double)(t4 - t3)));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD,
               "MUMPS factor memory estimate max-rank %.3f GiB, sum-ranks %.3f GiB, estimated factor entries %" PetscInt_FMT ", assembly peak max-rank %.3f GiB, sum-ranks %.3f GiB\n",
               (double)est_max_mb / 1024.0, (double)est_sum_mb / 1024.0, est_factor_entries,
               (double)(peak_max / (1024.0 * 1024.0 * 1024.0)), (double)(peak_sum / (1024.0 * 1024.0 * 1024.0))));
    PetscCall(MatDestroy(&factor));
    PetscCall(VecDestroy(&F));
    PetscCall(MatDestroy(&S));
    PetscCall(PetscFree(gids));
    PetscCall(PetscFree2(bdry, inner));
    PetscCall(PetscFree(shur_map));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  PetscCall(VecDuplicate(F, &U));
  PetscCall(VecDuplicate(F, &R));
  PetscCall(PetscTime(&t3));
  PetscCall(HmsCreateCoarseKSP(ctx, S, &coarse));
  PetscCall(PetscTime(&t4));
  PetscCall(VecSet(U, 0.0));
  PetscCall(KSPSolve(coarse, F, U));
  PetscCall(HmsCoarseComputeResidual(S, F, U, R, &rnorm, &fnorm));
  PetscCall(PetscTime(&t5));

  PetscCall(VecScatterCreateToAll(U, &scat, &Uall));
  PetscCall(VecScatterBegin(scat, U, Uall, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecScatterEnd(scat, U, Uall, INSERT_VALUES, SCATTER_FORWARD));
  PetscCall(VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, ms->global_p_nodes, &ufine));
  PetscCall(VecSet(ufine, 0.0));
  {
    const PetscScalar *ua = NULL;
    PetscScalar       *fa = NULL;
    PetscInt           lo, hi;
    PetscCall(VecGetArrayRead(Uall, &ua));
    PetscCall(VecGetOwnershipRange(ufine, &lo, &hi));
    PetscCall(VecGetArray(ufine, &fa));
    for (PetscInt gid = lo; gid < hi; gid++) if (shur_map[gid] >= 0) fa[gid - lo] = ua[shur_map[gid]];
    PetscCall(VecRestoreArray(ufine, &fa));

    for (PetscInt li = 0; li < local_cells; li++) {
      PetscInt      c = MsOwnedCellId(li, rank, size, ms->ms_spatial_cell_partition, cell_begin);
      PetscInt      cx = c % ms->coarse_x, tmp = c / ms->coarse_x, cz = tmp % ms->coarse_z, cy = tmp / ms->coarse_z;
      PetscScalar  *cb = NULL, *ui = NULL;
      PetscInt     *rows = NULL;
      PetscCall(BuildCellBoundaryGids(ms, shur_map, bdry, cx, cy, cz, gids));
      PetscCall(PetscMalloc2(ms->local_p_bdry, &cb, ms->local_p_inner, &rows));
      for (PetscInt q = 0; q < ms->local_p_bdry; q++) cb[q] = ua[gids[q]];
      PetscCall(FineSchurRecoverCell(ctx, ms, cx, cy, cz, cb, &ui, &recover_solve));
      for (PetscInt iy = 1; iy < ms->by; iy++) {
        for (PetscInt iz = 1; iz < ms->bz; iz++) {
          for (PetscInt ix = 1; ix < ms->bx; ix++) {
            PetscInt ii = LocalFineInteriorNodeId(ms, ix, iy, iz);
            rows[ii] = NodeId(ctx, cx * ms->bx + ix, cy * ms->by + iy, cz * ms->bz + iz);
          }
        }
      }
      PetscCall(VecSetValues(ufine, ms->local_p_inner, rows, ui, INSERT_VALUES));
      PetscCall(PetscFree(ui));
      PetscCall(PetscFree2(cb, rows));
    }
    PetscCall(VecRestoreArrayRead(Uall, &ua));
  }
  PetscCall(VecAssemblyBegin(ufine));
  PetscCall(VecAssemblyEnd(ufine));
  PetscCall(PetscTime(&t6));
  PetscCallMPI(MPI_Reduce(&recover_solve, &recover_solve_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));

  if (check_full) {
    PetscCall(AssembleSystem(ctx, &Afull, &bfull));
    PetscCall(VecDuplicate(bfull, &rfull));
    PetscCall(MatMult(Afull, ufine, rfull));
    PetscCall(VecAXPY(rfull, -1.0, bfull));
    PetscCall(VecNorm(rfull, NORM_2, &full_rnorm));
    PetscCall(VecNorm(bfull, NORM_2, &full_bnorm));
  }
  if (ctx->fine_reference_file[0]) {
    Vec       reference = NULL;
    PetscReal reference_norm = 0.0, difference_norm = 0.0;
    PetscCall(LoadVecBinary(PETSC_COMM_WORLD, ms->global_p_nodes, ctx->fine_reference_file, &reference));
    PetscCall(VecNorm(reference, NORM_2, &reference_norm));
    PetscCall(VecAXPY(reference, -1.0, ufine));
    PetscCall(VecNorm(reference, NORM_2, &difference_norm));
    reference_rel = reference_norm > 0.0 ? difference_norm / reference_norm : difference_norm;
    PetscCall(VecDestroy(&reference));
  }
  if (ctx->fine_solution_file[0]) PetscCall(SaveVecBinary(ufine, ctx->fine_solution_file));
  if (ctx->slice_prefix[0]) PetscCall(WriteSourceSlicesPrefix(ufine, ctx));
  PetscCall(PetscMemoryGetMaximumUsage(&peak));
  PetscCallMPI(MPI_Reduce(&peak, &peak_max, 1, MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD));
  PetscCallMPI(MPI_Reduce(&peak, &peak_sum, 1, MPI_DOUBLE, MPI_SUM, 0, PETSC_COMM_WORLD));

  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Fine exact-interface Schur direct solve\n"));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
             "fine mesh ne_y ne_z ne_x %" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT ", fine dofs %" PetscInt_FMT ", coarse cells %" PetscInt_FMT " %" PetscInt_FMT " %" PetscInt_FMT ", H/h %" PetscInt_FMT ", interface dofs %" PetscInt_FMT " (%.3f%%), local interior/boundary %" PetscInt_FMT "/%" PetscInt_FMT "\n",
             ctx->ne_y, ctx->ne_z, ctx->ne_x, ms->global_p_nodes, ms->coarse_y, ms->coarse_z, ms->coarse_x, ms->bx,
             ms->shur_nodes, 100.0 * (double)ms->shur_nodes / (double)ms->global_p_nodes, ms->local_p_inner, ms->local_p_bdry));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
             "time(s) interface pattern %.3f, local eliminate+Schur %.3f (PARDISO max-rank %.3f), interface MatAssembly %.3f, interface MUMPS factor %.3f, interface solve+check %.3f, recovery %.3f (PARDISO max-rank %.3f)\n",
             (double)(tp1 - tp0), (double)(t1 - t0), (double)local_solve_max, (double)(t2 - t1), (double)(t4 - t3), (double)(t5 - t4), (double)(t6 - t5), (double)recover_solve_max));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD,
             "interface |SU-F|/|F| %.6e, full |Au-b|/|b| %.6e, relative difference vs fine reference %.6e, peak memory max-rank %.3f GiB, sum-ranks %.3f GiB\n",
             (double)(fnorm > 0.0 ? rnorm / fnorm : rnorm), (double)(check_full && full_bnorm > 0.0 ? full_rnorm / full_bnorm : -1.0),
             (double)reference_rel,
             (double)(peak_max / (1024.0 * 1024.0 * 1024.0)), (double)(peak_sum / (1024.0 * 1024.0 * 1024.0))));

  PetscCall(VecDestroy(&rfull)); PetscCall(VecDestroy(&bfull)); PetscCall(MatDestroy(&Afull));
  PetscCall(VecDestroy(&ufine)); PetscCall(VecScatterDestroy(&scat)); PetscCall(VecDestroy(&Uall));
  PetscCall(KSPDestroy(&coarse)); PetscCall(VecDestroy(&R)); PetscCall(VecDestroy(&U)); PetscCall(VecDestroy(&F)); PetscCall(MatDestroy(&S));
  PetscCall(PetscFree(gids)); PetscCall(PetscFree2(bdry, inner)); PetscCall(PetscFree(shur_map));
  PetscFunctionReturn(PETSC_SUCCESS);
}
