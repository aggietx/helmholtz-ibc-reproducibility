include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

wavelet_outer_comparison: src/wavelet_ems_hs_3d_petsc.o
	${CLINKER} -o $@ $^ ${PETSC_KSP_LIB}

src/wavelet_ems_hs_3d_petsc.o: src/wavelet_ems_hs_3d_petsc.c src/hmsfem_q1_base.c src/fine_schur_direct.inc src/continuous_outer.inc
	${PCC} -o $@ -c ${PCC_FLAGS} ${CFLAGS} ${CCPPFLAGS} $<
