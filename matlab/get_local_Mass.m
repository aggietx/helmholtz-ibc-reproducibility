function [element_Mass]=get_local_Mass(n_node,hx,hy)
% %% p=n_node-1;
% % p=1;n_node=p+1;hx=2;hy=2;[local_Mass]=get_local_Masslag(n_node,hx,hy);
% % p=2;n_node=p+1;hx=2;hy=2;[local_Mass]=get_local_Masslag(n_node,hx,hy);
% n_node=3;

p=n_node-1;
gll=1;

if gll==1 %%% gll quadrature
[~,wx]=xwlgl(n_node);
[~,wy]=xwlgl(n_node);
favalue_x=speye(n_node);%first derivative GLL matrix % % saved in each column 
favalue_y=speye(n_node);
   
else %%%% gauss quadrature
 [x0]=gllquad(p,-1,1);%%%% inptepolation point
 [y0]=gllquad(p,-1,1);
[x,wx]=gaussquad(p,-1,1);
[y,wy]=gaussquad(p,-1,1);    
favalue_x=poly_quadvalue(x0,x);
favalue_y=poly_quadvalue(y0,y);
end

% if gll==1
% [~,wx]=xwlgl(n_node);
% [~,wy]=xwlgl(n_node);
% favalue_x=speye(n_node);%first derivative GLL matrix % % saved in each column 
% favalue_y=speye(n_node);
% else
% p=n_node-1;
% [x,wx]=gaussquad(p,-1,1);
% [y,wy]=gaussquad(p,-1,1);
% favalue_x= getlagquadvalue(x,p+1);
% favalue_y= getlagquadvalue(y,p+1);
% end

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