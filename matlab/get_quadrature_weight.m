function x = get_quadrature_weight(np,a,b) 
if np <= 1 
  es='The number of the quadrature nodes should be greater than 1';
  disp(es); x = []; w = [];
  return 
elseif np==2
    x=[-1;1];
elseif np==3
    x=[ -1;
     0;
     1];
elseif np==4
     x=[  -1.0000;
   -0.4472;
    0.4472;
    1.0000];
elseif np==5
     x=[  -1.0000;
   -0.6547;
    0.0000;
    0.6547;
    1.0000];
elseif np==6
   x=[    -1.0000;
   -0.7651;
   -0.2852;
    0.2852;
    0.7651;
    1.0000];
elseif np==7
    x=[   -1.0000;
   -0.8302;
   -0.4688;
         0;
    0.4688;
    0.8302;
    1.0000];
else %%% evaluate
    x = get_quadrature_point_weight(np,a,b);
end
   bma=(b-a)*.5;
   bpa=(b+a)*.5;
   x=bma*x+bpa;
   
% if np <= 1 
%   es='The number of the quadrature nodes should be greater than 1';
%   disp(es); x = []; w = [];
%   return 
%   
% elseif np==2
%     x=[-1;1];w=[1;1];
%     
% else
% x=zeros(np,1);
% w=zeros(np,1);
% n=np-1;
% x(1)=-1;x(np)=1;
% x1=jacobi_roots(n-1,1,1);
% x(2:n)=x1;
% coef=2/(n*np);
% w=coef./(pnleg(x,n).^2);
% end

% if nargin ==3 
%    bma=(b-a)*.5;
%    bpa=(b+a)*.5;
%    x=bma*x+bpa;
%    w=w*bma;
% end

return
