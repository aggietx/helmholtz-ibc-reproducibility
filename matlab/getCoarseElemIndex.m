function gi = getCoarseElemIndex(i, j, Nx, Ny, nx, ny)
%(i,j) is the sequence of the coarse elements, the output is all the nodes
%insides this element
nnx=nx*Nx+1;nny=ny*Ny+1;
all_element=reshape(1:nnx*nny,nnx,nny)';
gi=all_element((j-1)*ny+1:j*ny+1,(i-1)*nx+1:i*nx+1);gi=gi';
gi=gi(:);
% return
% ind = 1;
% for jj=1:ny+1
%   for ii=1:nx+1
%     gi(ind) = (ii +  (i-1)*nx) + (jj-1 + (j-1)*ny)*(Nx*nx + 1);
%     ind = ind+1;
%   end
% end