function [ix,iy]=get_ix_iy1(p,nx0,ny0)
% p=2;nex=2;ney=2;[ix,iy]=get_ix_iy1(p,nex,ney);
% % p is the degree of polynomial
% % total_element=nex*ney;
% localnod=getnod(p,nx0,ny0);
% localelement_node=localnod';% % nodes in each element(;% % dimension ((p+1)^2,nx*ny))
localelement_node=getnod1(p,nx0,ny0);

total_deg=(p+1)^2;
ix=repmat(localelement_node,total_deg,1);%ix=ix(:);

nov1=localelement_node(:);temp_y=repmat(nov1,1,total_deg);
iy=temp_y';iy=iy(:);iy=reshape(iy,total_deg^2,nx0*ny0);
