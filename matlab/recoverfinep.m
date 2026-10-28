function umsfine=recoverfinep(Ny,Nx,nx0,ny0,msbasis,uc,node_matrixc,p)
nx=Nx*nx0;ny=Ny*ny0;
umsfine=zeros((p*nx+1)*(p*ny+1),1);
    nnx=p*nx+1;nny=p*ny+1;
    all_dof=reshape(1:nnx*nny,nnx,nny)';
for j=1:Ny
  for i=1:Nx

      Neindex=(j-1)*Nx+i;
      localdof=getCoarsenodeselementinside(i,j,Nx,Ny,p*nx0,p*ny0,all_dof);
      localcoarsedof=node_matrixc( (j-1)*p+1:j*p+1,(i-1)*p+1:i*p+1);
      localcoarsedof=localcoarsedof.';
      localcoarsedof=localcoarsedof(:);
      umsfine(localdof)=msbasis(:,:,Neindex)*uc(localcoarsedof);

  end
end