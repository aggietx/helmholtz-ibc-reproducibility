function [ix,iy]=get_ix_iy(p,nov)
% p=2;nex=2;ney=2;[ix,iy]=get_ix_iy(p,nex,ney);
% % p is the degree of polynomial
% % total_element=nex*ney;
%%%% nov ( (p+1)^d,ne)

total_deg=(p+1)^2;
temp_x=repmat(nov,total_deg,1);ix=temp_x(:);
% whos temp_x
nov1=nov(:);temp_y=repmat(nov1,1,total_deg);
temp_y=temp_y';
iy=temp_y(:);
