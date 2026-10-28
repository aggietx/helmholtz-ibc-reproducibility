function [d] = derlgl(x,np) 
%DERLGL      Spectral (Legendre Gauss Lobatto) derivative matrix 
%
%    [d]=derlgl(x,np) returns the spectral Legendre Gauss Lobatto derivative
%    matrix d at the np LGL nodes x (on [-1,1]). np-1 is the polynomial
%    degree used. Legendre Gauss Lobatto (LGL) grid
%
%
% Input: x = array of LGL nodes on [-1,1] (computed by xwlgl)
%        np = number of LGL nodes (=n+1, n=polynomial interpolation degree)
%
% Output: D = spectral derivative matrix: formula  (2.3.28), pag. 80 CHQZ2
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$
 

d=zeros(np);
n=np-1; 
for j =1:np
   lnxj = pnleg(x(j),n);    
   for i = 1:np       
      if i ~= j          
         lnxi = pnleg(x(i),n);
         d(i,j) = lnxi/((x(i)-x(j))*lnxj);                 
      end       
   end 
end 
d(1,1) = -0.25*n*np; 
d(np,np) = -d (1,1); 
 
return
