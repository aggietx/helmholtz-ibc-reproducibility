function [x,y]=getcoor_sem(xl,xr,yl,yr,nx,ny,p)
%%% coordinates of all dof in a coarse element
%%%% for each fine-grid element I use its gll nodes
%%%% these nodes may not lies in the gll nodes of the coarse element


hx=(xr-xl)/nx;hy=(yr-yl)/ny;


x_1d=zeros(nx*p+1,1);
x_1d_initial=xl:hx:xr;
for i=1:nx
    a=x_1d_initial(i);b=x_1d_initial(i+1);
    quadrature_node = get_gll_weight_position(p,a,b); 
    x_1d((i-1)*p+1:i*p+1,1)=quadrature_node;
end
x=repmat(x_1d,ny*p+1,1);

y_1d=zeros(ny*p+1,1);
y_1d_initial=yl:hy:yr;
for i=1:ny
    a=y_1d_initial(i);b=y_1d_initial(i+1);
    quadrature_node = get_gll_weight_position(p,a,b); 
    y_1d((i-1)*p+1:i*p+1,1)=quadrature_node;
end
y_temp=repmat(y_1d,1,nx*p+1);y_temp=y_temp';y=y_temp(:);

