function imagescsquare(u)
u=real(u);
if size(u,1)==1
    u=u';
end
n=sqrt(size(u,1));
if size(u,1)==size(u,2)
    figure();imagesc(u)
else
    figure();imagesc(reshape(u,n,n)');colorbar
end
setgca(16);colormap jet;