function [element_Mass]=get_element_gauss_Mass(n_node,hx,hy)
%%%% p=1;n_node=p+1;hx=1;hy=1;[element_Mass]=get_element_gauss_Mass(n_node,hx,hy);
%%%%% legrendre polynomial with gauss mass
p=n_node-1;
[x,wx]=gaussquad(p,-1,1);
[y,wy]=gaussquad(p,-1,1);
favalue_x= getlen_gaussquadvalue(x,p+1);
favalue_y= getlen_gaussquadvalue(y,p+1);

jacx=hx/2;
jacy=hy/2;
jac=jacx*jacy;


npdx=length(wx);%%% equal n_node
npdy=length(wy);
mn=npdx*npdy;
A=zeros(mn,mn);
for iy1=1:npdy
     for ix1=1:npdx
          i1a=(iy1-1)*npdx+ix1; 
          for iy2=1:npdy
                 for ix2=1:npdx
                   i2a=(iy2-1)*npdx+ix2;  
                A(i1a,i2a)=   favalue_x(:,ix1)'*(diag(wx)*favalue_x(:,ix2))*...
                              favalue_y(:,iy1)'*(diag(wy)*favalue_y(:,iy2))*jac;
                 end
          end
      end
end

element_Mass=A;
% local_Mass=diag(A);