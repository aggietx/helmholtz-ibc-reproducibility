static PetscErrorCode W3FGMRES20(Mat A, Vec b, Vec x, W3PC *w, PetscInt cap,
                                PetscReal tol, PetscInt *its, PetscReal *rn,
                                KSPConvergedReason *reason)
{
  KSP ksp;
  PC pc;
  W3TrueResidualCtx check = {0};
  PetscReal bn;
  PetscInt restarts=0, cycle_its;
  PetscFunctionBeginUser;
  PetscCall(VecDuplicate(b,&check.solution));
  PetscCall(VecDuplicate(b,&check.residual));
  PetscCall(VecSet(check.solution,0.0));
  PetscCall(VecNorm(b,NORM_2,&bn));
  check.target_abs=tol*bn;
  PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
  PetscCall(KSPSetOperators(ksp,A,A));
  PetscCall(KSPSetType(ksp,KSPFGMRES));
  PetscCall(KSPGMRESSetRestart(ksp,20));
  PetscCall(KSPSetPCSide(ksp,PC_RIGHT));
  PetscCall(KSPSetNormType(ksp,KSP_NORM_UNPRECONDITIONED));
  PetscCall(KSPGetPC(ksp,&pc));
  PetscCall(PCSetType(pc,PCSHELL));
  PetscCall(PCShellSetContext(pc,w));
  PetscCall(PCShellSetApply(pc,W3PCApply));
  PetscCall(KSPSetConvergenceTest(ksp,W3TrueResidualConverged,&check,NULL));
  PetscCall(VecSet(x,0.0));
  *its=0; *rn=bn; *reason=KSP_DIVERGED_ITS;
  while (w->pc_calls<cap) {
    PetscInt remaining=cap-w->pc_calls;
    check.last_abs=PETSC_MAX_REAL; check.max_it=remaining;
    PetscCall(KSPSetInitialGuessNonzero(ksp,restarts ? PETSC_TRUE : PETSC_FALSE));
    PetscCall(KSPSetTolerances(ksp,tol,PETSC_DEFAULT,PETSC_DEFAULT,remaining));
    PetscCall(KSPSolve(ksp,b,x));
    PetscCall(KSPGetIterationNumber(ksp,&cycle_its));
    *its+=cycle_its;
    PetscCall(KSPGetConvergedReason(ksp,reason));
    PetscCall(MatMult(A,x,check.residual));
    PetscCall(VecAXPY(check.residual,-1.0,b));
    PetscCall(VecNorm(check.residual,NORM_2,rn));
    if (*rn<=tol*bn) { *reason=KSP_CONVERGED_RTOL; break; }
    if (*reason<0 || restarts>=8) break;
    restarts++;
  }
  PetscCall(KSPDestroy(&ksp));
  PetscCall(VecDestroy(&check.solution));
  PetscCall(VecDestroy(&check.residual));
  PetscFunctionReturn(PETSC_SUCCESS);
}
