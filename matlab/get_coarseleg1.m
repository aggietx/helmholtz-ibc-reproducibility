function polyc=get_coarseleg1(x,y,p)

ndof=size(x,1);
a = xwlgl0(p+1) ;%%% a is the gll points in (-1,1)

b=zeros(p+1,1);
for i=1:p+1
    f=1;
    for j=1:p+1 
        if j~=i
        f=f*(a(j)-a(i));
        end
    end
    b(i)=1/f;
end

F=zeros(ndof,(p+1)^2);
for j=1:p+1
    for i=1:p+1
        s=ones(ndof,1);
          for jj=1:p+1 
              if jj~=j 
          s=s.*(a(jj)-y);
              end
          end
              for ii=1:p+1
        if ii~=i
        s=s.*(a(ii)-x);
        end
              end
       s=b(i)*b(j)*s;
       F(:,(j-1)*(p+1)+i)=s;    
    end
end
polyc=F;

