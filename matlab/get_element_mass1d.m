function [local_Mass]=get_element_mass1d(n_node,hx)
%%%% n_node is the order+1
[~,wx]=get_gll_weight_position(n_node-1);

deri_x=speye(n_node);%first derivative GLL matrix % % saved in each column 


jac=hx/2;


npdx=length(wx);%%% equal n_node

mn=npdx;
ww=wx;
A=zeros(mn,mn);%
for ki=1:1
    inde=((ki-1)*npdx+1:ki*npdx);
    A(inde,inde)=deri_x'*(diag(ww(inde)))*jac;
end


local_Mass=A;
