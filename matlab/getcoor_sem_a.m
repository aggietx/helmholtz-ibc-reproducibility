function [x,y,x_1d,y_1d]=getcoor_sem_a(xl,xr,yl,yr,nx0,ny0,p)
%%%% coarse gll nodes then uniform mesh between the gll nodes
%%%% reference domain is [-1,1]*[-1,1]

quadrature_node = get_quadrature_weight(p+1,xl,xr);


x_1d=zeros(nx0*p+1,1);
x_1d(1:nx0:end)=quadrature_node;
for i=1:p
    xl=quadrature_node(i);
    xr=quadrature_node(i+1);
    
    x_1d( (i-1)*nx0+1:i*nx0+1)=xl:(xr-xl)/nx0:xr;
end
x=repmat(x_1d,ny0*p+1,1);

quadrature_node = get_quadrature_weight(p+1,yl,yr);
y_1d=zeros(ny0*p+1,1);
y_1d(1:ny0:end)=quadrature_node;
for i=1:p
    yl=quadrature_node(i);
    yr=quadrature_node(i+1);
    y_1d( (i-1)*ny0+1:i*ny0+1)=yl:(yr-yl)/ny0:yr;
end
y_temp=repmat(y_1d,1,nx0*p+1);y_temp=y_temp';y=y_temp(:);
