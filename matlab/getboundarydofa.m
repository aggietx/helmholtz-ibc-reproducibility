function [boundarynode,bd_nodes_pair]=getboundarydofa(nx,ny)
    
bd_nodes_pair=zeros(2*(nx+ny),2);
for i=1:nx
    bd_nodes_pair(i,:)=[i,i+1];
end
for i=nx+1:nx+ny
    bd_nodes_pair(i,:)=[(nx+1)*(i-nx),(nx+1)*(i+1-nx)];
end
for i=nx+ny+1:2*nx+ny
    bd_nodes_pair(i,:)=[(nx+1)*(ny+1)+nx+ny+1-i,(nx+1)*(ny+1)+nx+ny-i];
end
for i=2*nx+ny+1:2*(nx+ny)
    bd_nodes_pair(i,:)=[nx*(2*(ny+nx)+2)+nx+2*ny+2-(nx+1)*i,nx*(2*(nx+ny)+2)+nx+2*ny+2-(nx+1)*(i+1)];
end
boundarynode=bd_nodes_pair(:,1);