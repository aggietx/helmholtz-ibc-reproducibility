function localstiff1d=get_element_stiff1d_gauss(p,hx)
%%%% legendre poly with gauss quad
n_node=1+p;
[x,wx]=gaussquad(p,-1,1);
derix =derlen_gaussquadvalue(x,p+1); 

jacx=2/hx; 
localstiff1d=zeros(n_node,n_node);

for i1=1:n_node
    for i2=1:n_node
        localstiff1d(i1,i2)= derix(:,i1)'*(diag(wx)*derix(:,i2))*jacx;
    end
end