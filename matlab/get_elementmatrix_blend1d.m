function [element_Mass1d,elementstiff1d]=get_elementmatrix_blend1d(hx,tao,p)
element_Mass1dgll=get_element_mass1d(p+1,hx);
elementstiff1dgll=get_element_stiff1d(p,hx);
element_Mass1d_g=get_element_mass1d_gauss(p,hx);
elementstiff1d_g=get_element_stiff1d_gauss(p,hx);

element_Mass1d=tao*element_Mass1dgll+(1-tao)*element_Mass1d_g;
elementstiff1d=tao*elementstiff1dgll+(1-tao)*elementstiff1d_g;

