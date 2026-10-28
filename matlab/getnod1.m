function nod=getnod1(p,nx,ny)
%%%% the result is just the transpose of getnod
nod=zeros((p+1)^2,nx*ny);
all_nodes=1:(nx*p+1)*(ny*p+1);
node_matrix=reshape(all_nodes,nx*p+1,ny*p+1)';
for j=1:ny
    for i=1:nx
        temp_nodes=node_matrix( (j-1)*p+1:j*p+1,(i-1)*p+1:i*p+1);
        temp_nodes=temp_nodes';temp_nodes=temp_nodes(:);
        nod(:,(j-1)*nx+i)=temp_nodes;
    end
end