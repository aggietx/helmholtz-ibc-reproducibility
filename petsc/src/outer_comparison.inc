/* Standard right-preconditioned BiCGStab, accumulating the physical solution.
   True-residual checks require MatMult only, including the alpha half-step. */
static PetscErrorCode W3BiCGStab(Mat A, Vec b, Vec x, W3PC *w, PetscInt cap,
                               PetscReal tol, PetscInt *its, PetscReal *rn,
                               KSPConvergedReason *reason)
{
  Vec *v, r, shadow, p, ap, s, as, bp, bs, check;
  PC pc;
  PetscReal bn, asnorm;
  PetscLogDouble start, now;
  PetscScalar rho, oldrho = 1.0, alpha = 1.0, omega = 1.0, denom, beta;
  PetscFunctionBeginUser;
  PetscCall(VecDuplicateVecs(b, 9, &v));
  r=v[0]; shadow=v[1]; p=v[2]; ap=v[3]; s=v[4]; as=v[5]; bp=v[6]; bs=v[7]; check=v[8];
  PetscCall(PCCreate(PETSC_COMM_WORLD, &pc));
  PetscCall(PCSetOperators(pc, A, A));
  PetscCall(PCSetType(pc, PCSHELL));
  PetscCall(PCShellSetContext(pc, w));
  PetscCall(PCShellSetApply(pc, W3PCApply));
  PetscCall(PCSetUp(pc));
  PetscCall(VecNorm(b, NORM_2, &bn));
  PetscCall(PetscTime(&start));
  PetscCall(VecSet(x, 0.0));
  PetscCall(VecCopy(b, r));
  PetscCall(VecCopy(r, shadow));
  *rn=bn; *its=0; *reason=KSP_DIVERGED_ITS;
  if (bn == 0.0) *reason=KSP_CONVERGED_ATOL;
  while (*rn > tol*bn && w->pc_calls < cap) {
    /* PETSc VecDot(a,b) is b^H a for complex scalars. */
    PetscCall(VecDot(r, shadow, &rho));
    if (PetscAbsScalar(rho) == 0.0 || PetscIsInfOrNanScalar(rho)) {
      *reason=KSP_DIVERGED_BREAKDOWN; break;
    }
    if (*its == 0) PetscCall(VecCopy(r, p));
    else {
      beta=(rho/oldrho)*(alpha/omega);
      PetscCall(VecAXPY(p, -omega, ap));
      PetscCall(VecScale(p, beta));
      PetscCall(VecAXPY(p, 1.0, r));
    }
    PetscCall(PCApply(pc, p, bp));
    PetscCall(MatMult(A, bp, ap));
    PetscCall(VecDot(ap, shadow, &denom));
    if (PetscAbsScalar(denom) == 0.0 || PetscIsInfOrNanScalar(denom)) {
      *reason=KSP_DIVERGED_BREAKDOWN; break;
    }
    alpha=rho/denom;
    PetscCall(VecWAXPY(s, -alpha, ap, r));
    PetscCall(VecAXPY(x, alpha, bp));
    (*its)++;
    PetscCall(MatMult(A, x, check));
    PetscCall(VecAYPX(check, -1.0, b));
    PetscCall(VecNorm(check, NORM_2, rn));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "BICGSTAB PC=%" PetscInt_FMT " residual=%.12e\n", w->pc_calls, (double)(*rn/bn)));
    PetscCall(PetscTime(&now));
    if (*rn <= tol*bn) { *reason=KSP_CONVERGED_RTOL; break; }
    if (PetscIsInfOrNanReal(*rn) || *rn > 1e8*bn) { *reason=KSP_DIVERGED_DTOL; break; }
    if (w->pc_calls >= cap) break;
    if (now-start > 600.0) break;
    PetscCall(PCApply(pc, s, bs));
    PetscCall(MatMult(A, bs, as));
    PetscCall(VecNorm(as, NORM_2, &asnorm));
    PetscCall(VecDot(s, as, &denom));
    if (asnorm == 0.0 || PetscIsInfOrNanReal(asnorm)) { *reason=KSP_DIVERGED_BREAKDOWN; break; }
    omega=denom/(asnorm*asnorm);
    PetscCall(VecAXPY(x, omega, bs));
    PetscCall(VecWAXPY(r, -omega, as, s));
    PetscCall(MatMult(A, x, check));
    PetscCall(VecAYPX(check, -1.0, b));
    PetscCall(VecNorm(check, NORM_2, rn));
    if (*rn <= tol*bn) { *reason=KSP_CONVERGED_RTOL; break; }
    if (PetscIsInfOrNanReal(*rn) || *rn > 1e8*bn) { *reason=KSP_DIVERGED_DTOL; break; }
    if (PetscAbsScalar(omega) == 0.0 || PetscIsInfOrNanScalar(omega)) { *reason=KSP_DIVERGED_BREAKDOWN; break; }
    oldrho=rho;
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "BICGSTAB PC=%" PetscInt_FMT " residual=%.12e\n", w->pc_calls, (double)(*rn/bn)));
    PetscCall(PetscTime(&now));
    if (now-start > 600.0) break;
  }
  PetscCall(PCDestroy(&pc));
  PetscCall(VecDestroyVecs(9, &v));
  PetscFunctionReturn(PETSC_SUCCESS);
}
