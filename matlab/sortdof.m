function [allbddof,allinnerdof]=sortdof(nx,ny,p)
%%%% note 2d and 3d are different
alldof=reshape( 1: (nx*p+1)*(ny*p+1),nx*p+1,ny*p+1).';
allinnerdof=zeros((p-1)^2,nx*ny);
for j=1:ny
    for i=1:nx
        local_dof=alldof( (j-1)*p+1:j*p+1,(i-1)*p+1:i*p+1);
        localinnderdof=local_dof(2:end-1,2:end-1);
        localinnderdof=localinnderdof.';
       allinnerdof(:,(j-1)*nx+i)=localinnderdof(:);
    end
end

allinnerdof=allinnerdof(:);
allbddof=setdiff(1: (nx*p+1)*(ny*p+1),allinnerdof);
allbddof=allbddof.';