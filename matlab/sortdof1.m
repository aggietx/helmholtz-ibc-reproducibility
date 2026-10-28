function [allbddof,allinnerdof]=sortdof1(Nx,Ny,nx,ny,p)
%%%% for fine problem
alldof=reshape( 1: (Nx*nx*p+1)*(Ny*ny*p+1),Nx*nx*p+1,Ny*ny*p+1).';
allinnerdof=zeros((p*nx-1)*(p*ny-1),Nx*Ny);
for j=1:Ny
    for i=1:Nx
        local_dof=alldof( (j-1)*ny*p+1:j*ny*p+1,(i-1)*nx*p+1:i*nx*p+1);
        localinnderdof=local_dof(2:end-1,2:end-1);
        localinnderdof=localinnderdof.';
       allinnerdof(:,(j-1)*Nx+i)=localinnderdof(:);
    end
end

allinnerdof=allinnerdof(:);
allbddof=setdiff(1:(Nx*nx*p+1)*(Ny*ny*p+1),allinnerdof);
allbddof=allbddof.';