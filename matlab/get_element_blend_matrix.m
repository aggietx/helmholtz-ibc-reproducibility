
function [element_stiffness_matrix,element_stiffness_matrix_xx,element_stiffness_matrix_yy,...
    element_mass_matrix,element_Mass1d,element_stiff1d]=get_element_blend_matrix(hx,hy,p,tao)

[element_stiffness_matrix_g,element_stiffness_matrix_xx_g,element_stiffness_matrix_yy_g]=get_element_gauss_stiff(p+1,hx,hy);
[element_mass_matrix_g]=get_element_gauss_Mass(p+1,hx,hy);

[element_stiffness_matrix_gll,element_stiffness_matrix_xx_gll,element_stiffness_matrix_yy_gll]=get_elemet_stiff_matrix(p+1,hx,hy);
[element_mass_matrix_gll]=get_element_massmatrix(p+1,hx,hy);

element_Mass1d_g=get_element_mass1d_gauss(p,hx);
element_stiff1d_g=get_element_stiff1d_gauss(p,hx);
element_Mass1d_gll=get_element_mass1d(p+1,hx);%%% for 1d oscillating boundary condition
element_stiff1d_gll=get_element_stiff1d(p,hx);%%% for 1d oscillating boundary condition

element_stiffness_matrix_xx=tao*element_stiffness_matrix_xx_gll+(1-tao)*element_stiffness_matrix_xx_g;
element_stiffness_matrix_yy=tao*element_stiffness_matrix_yy_gll+(1-tao)*element_stiffness_matrix_yy_g;
element_mass_matrix=tao*element_mass_matrix_gll+(1-tao)*element_mass_matrix_g;
element_stiffness_matrix=tao*element_stiffness_matrix_gll+(1-tao)*element_stiffness_matrix_g;
element_Mass1d=tao*element_Mass1d_gll+(1-tao)*element_Mass1d_g;
element_stiff1d=tao*element_stiff1d_gll+(1-tao)*element_stiff1d_g;
