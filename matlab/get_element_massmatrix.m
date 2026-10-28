function [element_mass_matrix]=get_element_massmatrix(n_node,hx,hy)
%%% note this matrix is diagonal
% %% p=n_node-1;
% % n_node=3;hx=2;hy=2;[local_Mass]=get_local_Mass(n_node,hx,hy);
% %n_node=4;hx=2;hy=2;[local_Mass]=get_local_Mass(n_node,hx,hy);
% n_node=3;
order=n_node-1;
[~,wx]=get_gll_weight_position(order);
[~,wy]=get_gll_weight_position(order);
deri_x=speye(n_node);%first derivative GLL matrix % % saved in each column 
deri_y=speye(n_node);
jacx=hx/2;
jacy=hy/2;
jac=jacx*jacy;


npdx=length(wx);%%% equal n_node
npdy=length(wy);
% jacyx=jacy/jacx;
% jacxy=jacx/jacy;
mn=npdx*npdy;
[wx1,wy1]=meshgrid(wx,wy); ww=wx1.*wy1; ww=ww'; ww=ww(:); clear wx1 wy1;

A=zeros(mn,mn);%A=0;
% B=zeros(mn,mn);%B=0;
for ki=1:npdy
    inde=((ki-1)*npdx+1:ki*npdx);
    A(inde,inde)=deri_x'*(diag(ww(inde))*deri_y)*jac;
end


element_mass_matrix=A;
% local_Mass=diag(A);