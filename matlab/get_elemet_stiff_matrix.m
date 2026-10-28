function [elemet_stiff_xx_yy,elemet_stiff_xx,elemet_stiff_yy,ww]=get_elemet_stiff_matrix(n_node,hx,hy)
%%%%% sem element stiff matrix
%%%% n_node that elemet_stiff_xx,elemet_stiff_yy must be used in pml layer

%%% node is the order+1
%%% example
% %n_node=4; p=n_node-1;hx=1;hy=1;[local_Astiff]=get_elemet_stiff_matrix(p,hx,hy);


order=n_node-1;
[x,wx] = get_gll_weight_position(order);
[y,wy] = get_gll_weight_position(order);
deri_x=derlgl(x,n_node);%first derivative GLL matrix % % saved in each column 
deri_y=derlgl(y,n_node);
jacx=hx/2;
jacy=hy/2;
% % jac=jacx*jacy;

npdx=length(wx);
npdy=length(wy);
jacyx=jacy/jacx;
jacxy=jacx/jacy;
mn=npdx*npdy;
[wx1,wy1]=meshgrid(wx,wy); ww=wx1.*wy1; ww=ww'; ww=ww(:); clear wx1 wy1;

elemet_stiff_xx=zeros(mn,mn);%A=0;
elemet_stiff_yy=zeros(mn,mn);%B=0;
for ki=1:npdy
    inde=((ki-1)*npdx+1:ki*npdx);
    elemet_stiff_xx(inde,inde)=deri_x'*(diag(ww(inde))*deri_x)*jacyx;
end

for i=1:npdx
    inde=(i:npdx:mn);
    elemet_stiff_yy(inde,inde)=deri_y'*(diag(ww(inde))*deri_y)*jacxy;
end
elemet_stiff_xx_yy=elemet_stiff_xx+elemet_stiff_yy;
