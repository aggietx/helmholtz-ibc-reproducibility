function gi=getCoarsenodeselementinside(i,j,Nx,Ny,nx0,ny0,all_dof)
    
%     nnx=nx0*Nx+1;nny=ny0*Ny+1;
%     all_dof=reshape(1:nnx*nny,nnx,nny)';
    gi=all_dof((j-1)*ny0+1:j*ny0+1,(i-1)*nx0+1:i*nx0+1);

    gi=gi';
    gi=gi(:);gi=gi';
    
    
