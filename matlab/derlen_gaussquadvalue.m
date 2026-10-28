function [d,b] =derlen_gaussquadvalue(x,np) 
%%%%  legendre polynomial derivative matrix (each column stores a vector
%%%%% deri function's value in the quad points)
%%%%% x is gauss quad points on the reference interval [-1,1];
% [x,~]=gaussquad(p,-1,1);

p=np-1;
%% legendre polynomial with uniform dof distribution
a = get_gll_weight_position(p);
b=zeros(p+1,1);
for i=1:p+1
    f=1;
    for j=1:p+1 
        if j~=i
        f=f*(a(i)-a(j));
        end
    end
    b(i)=1/f;
end

%% deri values of lagrange polynomial on gauss points
d=zeros(np);
% p=np-1; 

for j =1:np
 
   s=zeros(np,1);
   for i = 1:np       
      if i ~= j   
          s1=ones(np,1);
         for jj=1:np
             if jj ~= i && jj~=j
             s1=-1*s1.*(a(jj)-x);
             end
         end
         s=s+s1;
      end 
      
   end 
   d(:,j)=s*b(j);

end 

%    d=d*(-1)^(np+1);
 
return
