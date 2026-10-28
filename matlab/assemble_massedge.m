function [Medge,Ml,Mr,Mb,Mt]=assemble_massedge(coeff,element_Mass1d,p)
ny=size(coeff,1);
nx=size(coeff,2);
dof=reshape(1: (p*ny+1)*(p*nx+1), (p*nx+1),(p*ny+1))';
ndof=(p*ny+1)*(p*nx+1);
element_Mass1d=element_Mass1d(:);
%% left
edgedof=dof(:,1);
nod=zeros(p+1,ny);
nod(1,:)=edgedof(1:p:end-1);
for i=2:p+1
nod(i,:)=edgedof(i:p:end);
end
ix=kron(nod,ones(p+1,1));iy=repmat(nod,p+1,1);
im=zeros( (p+1)^2,ny);
iix=1;
for iiy=1:ny
im(:,iiy)=coeff(iiy,iix)*element_Mass1d;
end

Ml=sparse(ix,iy,im,ndof,ndof);

%% right
edgedof=dof(:,end);
nod=zeros(p+1,ny);
nod(1,:)=edgedof(1:p:end-1);
for i=2:p+1
nod(i,:)=edgedof(i:p:end);
end
ix=kron(nod,ones(p+1,1));iy=repmat(nod,p+1,1);
im=zeros( (p+1)^2,ny);
iix=nx;
for iiy=1:ny
im(:,iiy)=coeff(iiy,iix)*element_Mass1d;
end

Mr=sparse(ix,iy,im,ndof,ndof);

%% bottom
edgedof=dof(1,:);
nod=zeros(p+1,nx);
nod(1,:)=edgedof(1:p:end-1);
for i=2:p+1
nod(i,:)=edgedof(i:p:end);
end
ix=kron(nod,ones(p+1,1));iy=repmat(nod,p+1,1);
im=zeros( (p+1)^2,nx);
iiy=1;
for iix=1:nx
im(:,iix)=coeff(iiy,iix)*element_Mass1d;
end

Mb=sparse(ix,iy,im,ndof,ndof);

%% bottom
edgedof=dof(end,:);
nod=zeros(p+1,nx);
nod(1,:)=edgedof(1:p:end-1);
for i=2:p+1
nod(i,:)=edgedof(i:p:end);
end
ix=kron(nod,ones(p+1,1));iy=repmat(nod,p+1,1);
im=zeros( (p+1)^2,nx);
iiy=ny;
for iix=1:nx
im(:,iix)=coeff(iiy,iix)*element_Mass1d;
end

Mt=sparse(ix,iy,im,ndof,ndof);
Medge=Mr+Mt+Ml+Mb;