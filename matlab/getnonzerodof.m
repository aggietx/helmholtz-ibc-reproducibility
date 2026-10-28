function [nonzerodof,zerodof]=getnonzerodof(nx_total,ny_total,k,nod)

%%% this code is to get nonzero and zero dof
% % % zeroindex= find(k.'~=0);
p=sqrt(size(nod,2))-1;
nonzeroindex= k.'~=0;
nonzerodof1=nod(nonzeroindex,:);
nonzerodof=unique(nonzerodof1);
zerodof=setdiff(1:(p*nx_total+1)*(p*ny_total+1),nonzerodof);
nonzerodof=sort(nonzerodof);
zerodof=sort(zerodof);
zerodof=zerodof.';