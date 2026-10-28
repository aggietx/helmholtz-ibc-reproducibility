function [p] = pnleg (x, n) 
%PNLEG      Evaluates Legendre polynomial of degree n at x
%
%    [p]=pnleg(x,n) returns the evaluation of the Legendre
%    polynomial of degree n at x, using the three
%    term relation.
%    (formula (2.3.3), pag. 75 CHQZ2)
%
% Input: x = scalar or column array
%        n = degree of Legendre polynomial to be evalueted
%
% Output: p = scalar or column array (same dimension as x)
%             containing the  evaluation of L_n at x
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


m=length(x);
if m==1

p=0;
if n==0
p=1;
else
p1=1;  p2=x;  p3=p2; 
for k =1:n-1    
   p3=((2*k+1)*x*p2-k*p1)/(k+1); 
   p1=p2; 
   p2=p3; 
end 
p=p3; 
end

else
x=x(:);
p=zeros(m,1);
if n==0
p=ones(m,1);
else
p1=ones(m,1);  p2=x;  p3=p2;
for k =1:n-1
   p3=((2*k+1)*x.*p2-k*p1)/(k+1);
   p1=p2;
   p2=p3;
end
p=p3;
end

end 
return 
