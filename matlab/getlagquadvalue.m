function [favale,polycoeff] = getlagquadvalue(x,np) 
%%% get the quad value of lagrange functions (each column stores a vector
%%%%%function's value in the quad points)
%%%%% x is gauss quad points on the reference interval [-1,1];
% p=2;[x,~]=gaussquad(p,-1,1);[d,b] = getlagquadvalue(x,1+p);
% 

p=np-1;
%% lagrange polynomial with uniform dof distribution
a=-1:2/p:1;
polycoeff=zeros(p+1,1);
for i=1:p+1
    f=1;
    for j=1:p+1 
        if j~=i
        f=f*(a(i)-a(j));
        end
    end
    polycoeff(i)=1/f;
end

%% deri values of lagrange polynomial on gauss points
favale=zeros(np);
% p=np-1; 

for j =1:np
 
        s=ones(np,1);
          for jj=1:p+1 
              if jj~=j 
          s=s.*(a(jj)-x);
              end
          end
        favale(:,j)=polycoeff(j)*s;    

end 

% favale=favale*(-1)^(np+1);
 
return
