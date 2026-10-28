function z = dd_pre_hel_mibc_2level_hybrid_batched(phiA, Phi, ...
    Uc, Lc, Pc, Qc, dd, inc)
%DD_PRE_HEL_MIBC_2LEVEL_HYBRID_BATCHED Apply Q0(I-A*BRAS)+BRAS.

z1 = dd_pre_hel_mibc_ras_pou_batched(dd, inc);
coarse_solve = @(r0) Qc * (Uc \ (Lc \ (Pc * r0)));
z = z1 - Phi * coarse_solve(phiA * z1) + ...
    Phi * coarse_solve(Phi' * inc);
end
