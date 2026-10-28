function meshsquare(u)
u=real(u);
if size(u,2)==size(u,1)
figure();mesh(u)
else
if size(u,1)==1
    u=u';
end
n=sqrt(size(u,1));
figure();mesh(reshape(u,n,n)')
end