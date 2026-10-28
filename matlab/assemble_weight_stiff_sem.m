function [Astiff]=assemble_weight_stiff_sem(local_Astiff,element_model,ix,iy)
%%% ix iy from [ix,iy]=get_ix_iy(p,nex,ney)
%%% eg [ix,iy]=get_ix_iy(1,2,2);
%%%% element_model=[0.5,1;1.5,2];
% total_deg=size(local_Astiff,2)*size(local_Astiff,1);
% ngrid=size(element_model,1)*size(element_model,2);
element_model=element_model.';
coeff=element_model(:);
% coeff_matrix=repmat(coeff.',total_deg,1);

% whos coeff_matrix


la=local_Astiff(:);
% lk_repeat=la(:,ones(ngrid,1));
% whos lk_repeat;

% ivalue=coeff_matrix.*lk_repeat;
ivalue=la*coeff.';
ivalue=ivalue(:);
Astiff=sparse(ix,iy,ivalue);