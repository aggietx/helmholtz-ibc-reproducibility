function [x,w] = xwlgl(np,a,b) 



if np <= 1 
  es='The number of the quadrature nodes should be greater than 1';
  disp(es); x = []; w = [];
  return 
  
elseif np==2
    x=[-1;1];w=[1;1];
    
else
x=zeros(np,1);
w=zeros(np,1);
n=np-1;
x(1)=-1;x(np)=1;
x1=jacobi_roots(n-1,1,1);
x(2:n)=x1;
coef=2/(n*np);
w=coef./(pnleg(x,n).^2);
end

if nargin ==3 
   bma=(b-a)*.5;
   bpa=(b+a)*.5;
   x=bma*x+bpa;
   w=w*bma;
end

return
