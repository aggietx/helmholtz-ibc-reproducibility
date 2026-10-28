function [x,flag] = jacobi_roots(n,alpha,beta) 
% JACOBI_ROOTS: Computes the n zeros of the Jacoby polynomial P_n^{(\alpha,\beta)}(x)
%   by Newton method and deflation process.
%
%    [x] = jacobi_roots(n,alpha,beta) 
%
% Input: n = polynomial degree
%        alpha,beta = parameters of Jacobi polynomial
%
% Output: x = zeros of Jacobi polynomial (column array of size n)
%         flag  = -1: The polynomial degree should be greater than 0
%
% Reference: CHQZ2 = C. Canuto, M.Y. Hussaini, A. Quarteroni, T.A. Zang,
%                    "Spectral Methods. Fundamentals in Single Domains"
%                    Springer Verlag, Berlin Heidelberg New York, 2006.

%   Written by Paola Gervasio
%   $Date: 2007/04/01$


flag=0;
if n < 1 
  es='The polynomial degree should be greater than 0';
  disp(es); flag = -1; x = []; 
return 
 
else
x=zeros(n,1);

x0=cos(pi/(2*n));
tol=1.e-14; kmax=15;
for j=1:n
 diff=tol+1;kiter=0;
    while kiter <=kmax & diff>=tol
        [p,pd]=jacobi_eval(x0,n,alpha,beta);
% deflation process q(x)=p(x)*(x-x_1)*... (x-x_{j-1})
% q(x)/q'(x)=p(x)/[p'(x)-p(x)*\sum_{i<j} 1/(x-x_i)]
        ss=sum(1./(x0-x(1:j-1)));
        x1=x0-p/(pd-ss*p);
        diff=abs(x1-x0);
        kiter=kiter+1;
        x0=x1;
    end
    x0=5.d-1*(x1+cos((2*(j+1)-1)*pi/(2*n)));
    x(j)=x1;
end 
x=sort(x);
end


return
