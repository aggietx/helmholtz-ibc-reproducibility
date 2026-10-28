function [p,pd] = jacobi_eval(x,n,alpha,beta) 
% JACOBI_EVAL: Evaluates Jacobi polynomial P_n^{(\alpha,\beta)} at x\in(-1,1)
%
%              (formula (2.5.3) pag. 92, CHQZ2)
%  [p,pd] = jacobi_eval(x,n,alpha,beta) 
%
% Input: x = scalar or one-dimensional array of length (m) 
%        n =  polynomial degree
%        alpha, beta= parameters of Jacoby polynomial
%
% Output: p(m,3) = [P_n^{(\alpha,\beta)}(x),
%                   P_(n-1)^{(\alpha,\beta)}(x),
%                   P_(n-2)^{(\alpha,\beta)}(x)];
%
%         pd(m,3) = [(P_n^{(\alpha,\beta)})'(x),
%                    (P_(n-1)^{(\alpha,\beta)})'(x),
%                    (P_(n-2)^{(\alpha,\beta)})'(x)];
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$

apb=alpha+beta; ab2=alpha^2-beta^2; 

nn=length(x);
if nn==1
% x is a scalar
p=1;   pd=0; 
p1=0;  pd1=0;
p2=0;  pd2=0;
if n == 0    
   return 
elseif n==1
p1 = p; p2=p1;
p = (alpha-beta+(apb+2)*x)*0.5; 

pd1 = pd; pd2=pd1;
pd = 0.5*(apb+2); 
else
p1 = p; p2=p1;
p = (alpha-beta+(apb+2)*x)*0.5; 

pd1 = pd; pd2=pd1;
pd = 0.5*(apb+2); 
for k = 1:n-1 
   k1=k+1; k2=k*2; k2ab=k2+alpha+beta;
   k2ab1=k2ab+1; k2ab2=k2ab1+1;
   p2=p1; p1=p; 
   pd2=pd1; pd1=pd; 
   a1 = 2*k1*(k1+apb)*k2ab; 
% Gamma(x+n)/Gamma(x) = (x+n-1) * ... * (x+1) * x 
   a21 = k2ab1*ab2;
   a22 = k2ab2*k2ab1*k2ab; 
   a3=2*(k+alpha)*(k+beta)*k2ab2;
   p = ((a21+a22*x)*p1-a3*p2)/a1; 
   pd= (a22*p1+(a21+a22*x)*pd1-a3*pd2)/a1; 
end 
end

else
% x is an array

[m1,m2]=size(x);
if m1<m2
    x=x';
end
m=max(m1,m2);
p=[ones(m,1),zeros(m,1),zeros(m,1)];   pd=zeros(m,3); 
if n == 0    
   return 
elseif n==1
p(:,2) = p(:,1); p(:,3)=p(:,2);
p(:,1) = (alpha-beta+(apb+2)*x)*0.5; 

pd(:,2) = pd(:,1); pd(:,3)=pd(:,2);
pd(:,1) = 0.5*(apb+2); 
else
p(:,2) = p(:,1); p(:,3)=p(:,2);
p(:,1) = (alpha-beta+(apb+2)*x)*0.5; 

pd(:,2) = pd(:,1); pd(:,3)=pd(:,2);
pd(:,1) = 0.5*(apb+2);     
for k = 1:n-1 
    k2=k*2; k2ab=k2+alpha+beta;
   p(:,3)=p(:,2); p(:,2)=p(:,1); 
   pd(:,3)=pd(:,2); pd(:,2)=pd(:,1); 
   a1 = 2*(k+1)*(k+apb+1)*k2ab; 
% Gamma(x+n)/Gamma(x) = (x+n-1) * ... * (x+1) * x 
   a21 = (k2ab+1)*ab2;
   a22 = (k2ab+2)*(k2ab+1)*k2ab; 
   a3=2*(k+alpha)*(k+beta)*(k2ab+2);
   p(:,1) = ((a21+a22*x).*p(:,2)-a3*p(:,3))/a1; 
   pd(:,1)= (a22*p(:,2)+(a21+a22*x).*pd(:,2)-a3*pd(:,3))/a1; 
end 
end

end
return 
 
 
