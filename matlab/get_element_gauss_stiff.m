function [A,Axx,Ayy]=get_element_gauss_stiff(n_node,hx,hy)
%%%% p=1;n_node=p+1;hx=1;hy=1;[A,Axx,Ayy]=get_element_gauss_stiff(n_node,hx,hy);
%%%%% legrendre polynomial with gauss mass
p=n_node-1;
[x,wx]=gaussquad(p,-1,1);
[y,wy]=gaussquad(p,-1,1);
favaluex= getlen_gaussquadvalue(x,p+1);
favaluey= getlen_gaussquadvalue(y,p+1);
derix =derlen_gaussquadvalue(x,p+1); 
deriy =derlen_gaussquadvalue(x,p+1);

jacx=hx/2;
jacy=hy/2;



npdx=length(wx);%%% equal n_node
npdy=length(wy);
mn=npdx*npdy;
Axx=zeros(mn,mn);
Ayy=zeros(mn,mn);
for iy1=1:npdy
     for ix1=1:npdx
          i1a=(iy1-1)*npdx+ix1; 
          for iy2=1:npdy
                 for ix2=1:npdx
                   i2a=(iy2-1)*npdx+ix2;  
                Axx(i1a,i2a)=   derix(:,ix1)'*(diag(wx)*derix(:,ix2))*...
                              favaluey(:,iy1)'*(diag(wy)*favaluey(:,iy2))*jacy/jacx;
                Ayy(i1a,i2a)=   favaluex(:,ix1)'*(diag(wx)*favaluex(:,ix2))*...
                              deriy(:,iy1)'*(diag(wy)*deriy(:,iy2))*jacx/jacy;                                       
                 end
          end
      end
end

A=Axx+Ayy;