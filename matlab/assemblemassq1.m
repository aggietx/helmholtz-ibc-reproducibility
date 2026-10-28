function [Mass_weight]=assemblemassq1(local_Mass,element_model,ix,iy)
% % local_Mass is a vector
%%% ix iy from [ix,iy]=get_ix_iy(p,nex,ney)
%%% eg [ix,iy]=get_ix_iy(1,2,2);
%%%% element_model=[0.5,1;1.5,2];
% total_deg=size(local_Mass,2)*size(local_Mass,1);
% ngrid=size(element_model,1)*size(element_model,2);
element_model=element_model.';
coeff=element_model(:);
% coeff_matrix=repmat(coeff',total_deg,1);


la=local_Mass(:);
% lk_repeat=la(:,ones(ngrid,1));
% ivalue=coeff_matrix.*lk_repeat;
ivalue=la*coeff.';ivalue=ivalue(:);
Mass_weight=sparse(ix,iy,ivalue); % % matrix  but is diagnoal
