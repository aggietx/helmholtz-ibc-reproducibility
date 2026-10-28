function pu1dx=generatelinearpu1d(nx,Nx,overe)
nnx=nx*Nx+1;
pu1dx=zeros(nnx,Nx);

pu1dx(1:nx+1-overe,1)=1;pu1dx(nx+1-overe:nx+1+overe,1)=1:-1/(2*overe):0;

for i=2:Nx-1
pu1dx((i-1)*nx+overe+1: i*nx+1-overe,i)=1;
pu1dx(i*nx+1-overe:i*nx+1+overe,i)=1:-1/(2*overe):0;
pu1dx((i-1)*nx+1-overe:(i-1)*nx+1+overe,i)=0:1/(2*overe):1;
end
pu1dx( (Nx-1)*nx+1+overe:end,Nx)=1;pu1dx((Nx-1)*nx+1-overe:(Nx-1)*nx+1+overe,Nx)=0:1/(2*overe):1;
