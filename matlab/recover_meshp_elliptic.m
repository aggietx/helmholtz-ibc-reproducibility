function umsfine=recover_meshp_elliptic(Ny,Nx,nx0,ny0,trailbasisall,uc,node_matrixc,p)
%%% only recover partial p
nx=Nx*nx0;ny=Ny*ny0;
umsfine=zeros((1*nx+1)*(1*ny+1),1);
    nnx=1*nx+1;nny=1*ny+1;
    all_dof=reshape(1:nnx*nny,nnx,nny)';
for j=1:Ny
  for i=1:Nx

      Neindex=(j-1)*Nx+i;
      localdof=getCoarsenodeselementinside(i,j,Nx,Ny,1*nx0,1*ny0,all_dof);
% %       localdof=[localdof,localdof+(1*nx+1)*(1*ny+1)];
      localcoarsedof=node_matrixc( (j-1)*p+1:j*p+1,(i-1)*p+1:i*p+1);
      localcoarsedof=localcoarsedof.';
      localcoarsedof=localcoarsedof(:);
% % %       localcoarsedof=[localcoarsedof;localcoarsedof+(p*Nx+1)*(p*Ny+1)];
      umsfine(localdof)=trailbasisall(:,:,Neindex)*uc(localcoarsedof);

  end
end