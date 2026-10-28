function localstiff1d=get_element_mass1d_gauss(p,hx)
%%%% legendre poly with gauss quad
n_node=1+p;
[x,wx]=gaussquad(p,-1,1);
favaluex= getlen_gaussquadvalue(x,p+1);
% % derix =derlen_gaussquadvalue(x,p+1); 

jacx=hx/2; 
localstiff1d=zeros(n_node,n_node);

for i1=1:n_node
    for i2=1:n_node
        localstiff1d(i1,i2)= favaluex(:,i1)'*(diag(wx)*favaluex(:,i2))*jacx;
    end
end