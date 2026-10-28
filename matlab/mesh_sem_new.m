function [Coordinates,x,y,nod,n_global_nodes,nem,free_nodes,bdnodes,node_matrix,...
    bdbasisdof,innerbasisdof,nodinnerbasis]=mesh_sem_new(Lx,Ly,nx,ny,p)
%Supply the mesh information
nem=nx*ny;%The number of elements
all_nodes=1:(nx*p+1)*(ny*p+1);
node_matrix=reshape(all_nodes,nx*p+1,ny*p+1)';
nod=zeros(nem,(p+1)^2);
nodinnerbasis=zeros( (p-1)^2,nx*ny);
for j=1:ny
    for i=1:nx
        temp_nodes=node_matrix( (j-1)*p+1:j*p+1,(i-1)*p+1:i*p+1);
        temp_nodes=temp_nodes';temp_nodes=temp_nodes(:);
        nod((j-1)*nx+i,:)=temp_nodes;
        if p>1
            temp_innerdof=node_matrix( (j-1)*p+2:j*p,(i-1)*p+2:i*p);
            temp_innerdof=temp_innerdof.';temp_innerdof=temp_innerdof(:);
            ne=(j-1)*nx+i;
            nodinnerbasis(:,ne)=temp_innerdof;
        end
    end
end
innerbasisdof=sort(nodinnerbasis(:));
bdbasisdof=setdiff(all_nodes,innerbasisdof);


free_nodes=node_matrix(2:end-1,2:end-1);free_nodes=sort(free_nodes(:));
bdnodes=setdiff(all_nodes',free_nodes);

hx=Lx/nx;hy=Ly/ny;
n_global_nodes=(nx*p+1)*(ny*p+1);%The number of nodes

x_1d=zeros(nx*p+1,1);
x_1d_initial=0:hx:Lx;
for i=1:nx
    a=x_1d_initial(i);b=x_1d_initial(i+1);
    [quadrature_node]=get_quadrature_point_weight(p+1,a,b);
    x_1d((i-1)*p+1:i*p+1,1)=quadrature_node;
end
x=repmat(x_1d,ny*p+1,1);

y_1d=zeros(ny*p+1,1);
y_1d_initial=0:hy:Ly;
for i=1:ny
    a=y_1d_initial(i);b=y_1d_initial(i+1);
    [quadrature_node]=get_quadrature_point_weight(p+1,a,b);
    y_1d((i-1)*p+1:i*p+1,1)=quadrature_node;
end
y_temp=repmat(y_1d,1,nx*p+1);y_temp=y_temp';y=y_temp(:);

Coordinates=[x,y]; 
