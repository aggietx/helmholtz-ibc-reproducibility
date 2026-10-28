function [fvalue,polycoeff] =poly_quadvalue(a,x) 
%%%% polynomial value matrix (each column stores a vector
%%%%%  function's value in the quad points)
% example of x;[x,~]=gaussquad(3,-1,1);%% quad points
% % example of a; a = gllquad(2,-1,1);%%%% zero points of the polynomial,if a is uniform distributed then it is lagrange polynomial
%%% a determines polynomial, x determines polynomial value
%%%% size of a and x may be different

%% polynomial coefficient
p1=size(a,1)-1;
polycoeff=zeros(p1+1,1);
for i=1:p1+1
    f=1;
    for j=1:p1+1 
        if j~=i
        f=f*(a(j)-a(i));
        end
    end
    polycoeff(i)=1/f;
end

%% values of polynomial on  points to be evaluated(x)
p2=size(x,1)-1;
fvalue=zeros(p2+1,p1+1);
for j =1:p1+1
        s=ones(p2+1,1);
          for jj=1:p1+1 
              if jj~=j
          s=s.*(a(jj)-x);
              end
          end
        fvalue(:,j)=polycoeff(j)*s;  
end 

%    d=d*(-1)^(np+1);
 
return
