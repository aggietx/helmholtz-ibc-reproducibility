function A=assemblymatrix(ix,iy,elementmodel,local_matrix)
%%% ix and iy are matrix
elementmodel=elementmodel.';
coeff=elementmodel(:);



la=local_matrix(:);

ivalue=la*coeff.';
A=sparse(ix,iy,ivalue);