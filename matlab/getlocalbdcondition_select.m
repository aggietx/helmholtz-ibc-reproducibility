function [localbd,localbd1]=getlocalbdcondition_select(nx,ny,type)
localbd=zeros(2*nx+2*ny,64);

% if mod(nx,16)~=0 || mod(ny,16)~=0
% disp('error, localnx and localny must be divided by 16')
% end
% bottom=1:nx+1;right=nx+1:nx+ny+1;
% top=nx+ny+1:2*nx+ny+1;left=[2*nx+ny+1:2*nx+2*ny,1];
%type=1;
% i=10; mesh(reshape(localbd1(:,i),nx+1,ny+1)')
if type==1 %%%% continuous boundary condition
    
localbd(1:nx+1,1)=1:-1/nx:0;localbd(2*nx+ny+1:end,1)=0:1/ny:1-1/ny;
localbd(1:nx+1,2)=0:1/nx:1;localbd(nx+1:nx+ny+1,2)=1:-1/ny:0;
localbd(nx+1:nx+ny+1,3)=0:1/ny:1;localbd(nx+ny+1:2*nx+ny+1,3)=1:-1/nx:0;
localbd(nx+ny+1:2*nx+ny+1,4)=0:1/nx:1;localbd(2*nx+ny+1:end,4)=1:-1/ny:1/ny;
    
a=0; n1=floor(nx/2);n2=nx-n1; localbd(a+1:a+n1+1,5)=0:1/n1:1;localbd(a+n1+1:a+nx+1,5)=1:-1/(n2):0;
    
a=nx;n1=floor(ny/2);n2=ny-n1; localbd(a+1:a+n1+1,6)=0:1/n1:1;localbd(a+n1+1:a+ny+1,6)=1:-1/n2:0;
   
a=nx+ny;n1=floor(nx/2);n2=nx-n1;localbd(a+1:a+n1+1,7)=0:1/n1:1;localbd(a+n1+1:a+nx+1,7)=1:-1/n2:0;
a=nx*2+ny;n1=floor(ny/2);n2=ny-n1;localbd(a+1:a+n1+1,8)=0:1/n1:1;localbd(a+n1+1:a+ny,8)=1:-1/n2:1/n2;
     
a=0;nt=floor(nx/2);n1=floor(nx/4);n2=nt-n1;localbd(a+1:a+n1+1,9)=0:1/n1:1;localbd(a+n1+1:a+nt+1,9)=1:-1/n2:0;
   
a=floor(nx/2);nt=nx-floor(nx/2);n1=floor(nx/4);n2=nt-n1;localbd(a+1:a+n1+1,10)=0:1/(n1):1;localbd(a+n1+1:a+nt+1,10)=1:-1/(n2):0;
 
a=nx;nt=floor(ny/2);n1=floor(ny/4);n2=nt-n1;localbd(a+1:a+n1+1,11)=0:1/(n1):1;localbd(a+n1+1:a+nt+1,11)=1:-1/(n2):0;

a=nx+floor(ny/2);nt=nx-floor(nx/2);n1=floor(nx/4);n2=nt-n1;localbd(a+1:a+n1+1,12)=0:1/(n1):1;localbd(a+n1+1:a+nt+1,12)=1:-1/(n2):0;

a=nx+ny;nt=floor(nx/2);n1=floor(nx/4);n2=nt-n1;localbd(a+1:a+n1+1,13)=0:1/(n1):1;localbd(a+n1+1:a+nt+1,13)=1:-1/(n2):0;
a=nx+ny+floor(nx/2);nt=floor(nx/2);n1=floor(nx/4);n2=nt-n1;localbd(a+1:a+n1+1,14)=0:1/(n1):1;localbd(a+n1+1:a+nt+1,14)=1:-1/(n2):0;

a=nx+ny+nx;nt=floor(ny/2);n1=floor(ny/4);n2=nt-n1;localbd(a+1:a+n1+1,15)=0:1/(n1):1;localbd(a+n1+1:a+nt+1,15)=1:-1/(n2):0;
a=nx+ny+nx+floor(ny/2);nt=floor(ny/2);n1=floor(ny/4);n2=nt-n1;localbd(a+1:a+n1+1,16)=0:1/(n1):1;localbd(a+n1+1:a+nt,16)=1:-1/(n2):1/(n2);

   if 1
a=0;localbd(a+1:a+nx/8+1,17)=0:1/(nx/8):1;localbd(a+nx/8+1:a+nx/4+1,17)=1:-1/(nx/8):0;
a=nx/4;localbd(a+1:a+nx/8+1,18)=0:1/(nx/8):1;localbd(a+nx/8+1:a+nx/4+1,18)=1:-1/(nx/8):0;
a=2*nx/4;localbd(a+1:a+nx/8+1,19)=0:1/(nx/8):1;localbd(a+nx/8+1:a+nx/4+1,19)=1:-1/(nx/8):0;
a=3*nx/4;localbd(a+1:a+nx/8+1,20)=0:1/(nx/8):1;localbd(a+nx/8+1:a+nx/4+1,20)=1:-1/(nx/8):0;

d=nx/16;
for i=1:8
a=(i-1)*2*d;id=32+i;localbd(a+1:a+d+1,id)=0:1/d:1;localbd(a+d+1:a+2*d+1,id)=1:-1/d:0;
end

a=nx;localbd(a+1:a+ny/8+1,21)=0:1/(ny/8):1;localbd(a+ny/8+1:a+ny/4+1,21)=1:-1/(ny/8):0;
a=nx+ny/4;localbd(a+1:a+ny/8+1,22)=0:1/(ny/8):1;localbd(a+ny/8+1:a+ny/4+1,22)=1:-1/(ny/8):0;
a=nx+ny/2;localbd(a+1:a+ny/8+1,23)=0:1/(ny/8):1;localbd(a+ny/8+1:a+ny/4+1,23)=1:-1/(ny/8):0;
a=nx+3*ny/4;localbd(a+1:a+ny/8+1,24)=0:1/(ny/8):1;localbd(a+ny/8+1:a+ny/4+1,24)=1:-1/(ny/8):0;
d=ny/16;
for i=1:8
a=nx+(i-1)*2*d;id=40+i;localbd(a+1:a+d+1,id)=0:1/d:1;localbd(a+d+1:a+2*d+1,id)=1:-1/d:0;
end

a=nx+ny;localbd(a+1:a+nx/8+1,25)=0:1/(nx/8):1;localbd(a+nx/8+1:a+nx/4+1,25)=1:-1/(nx/8):0;
a=nx+ny+nx/4;localbd(a+1:a+nx/8+1,26)=0:1/(nx/8):1;localbd(a+nx/8+1:a+nx/4+1,26)=1:-1/(nx/8):0;
a=nx+ny+2*nx/4;localbd(a+1:a+nx/8+1,27)=0:1/(nx/8):1;localbd(a+nx/8+1:a+nx/4+1,27)=1:-1/(nx/8):0;
a=nx+ny+3*nx/4;localbd(a+1:a+nx/8+1,28)=0:1/(nx/8):1;localbd(a+nx/8+1:a+nx/4+1,28)=1:-1/(nx/8):0;

d=nx/16;
for i=1:8
a=nx+ny+(i-1)*2*d;id=48+i;localbd(a+1:a+d+1,id)=0:1/d:1;localbd(a+d+1:a+2*d+1,id)=1:-1/d:0;
end
a=nx+ny+nx;localbd(a+1:a+ny/8+1,29)=0:1/(ny/8):1;localbd(a+ny/8+1:a+ny/4+1,29)=1:-1/(ny/8):0;
a=nx+ny+nx+ny/4;localbd(a+1:a+ny/8+1,30)=0:1/(ny/8):1;localbd(a+ny/8+1:a+ny/4+1,30)=1:-1/(ny/8):0;
a=nx+ny+nx+ny/2;localbd(a+1:a+ny/8+1,31)=0:1/(ny/8):1;localbd(a+ny/8+1:a+ny/4+1,31)=1:-1/(ny/8):0;
a=nx+ny+nx+3*ny/4;localbd(a+1:a+ny/8+1,32)=0:1/(ny/8):1;localbd(a+ny/8+1:a+ny/4,32)=1:-1/(ny/8):1/(ny/8);

d=ny/16;
for i=1:8
a=2*nx+ny+(i-1)*2*d;id=56+i;localbd(a+1:a+d+1,id)=0:1/d:1;
localbd(a+d+1:a+2*d,id)=1:-1/d:1/d;
end
    end
else
    localbd(1:nx,1)=1;
    localbd(nx+1:nx+ny,2)=1;
    localbd(nx+ny+1:2*nx+ny,3)=1;
    localbd(2*nx+ny+1:end,4)=1;
    localbd(1:nx/2,5)=1;localbd(nx/2+1:nx,5)=-1;
    localbd(1+nx:ny/2+nx,6)=1;localbd(ny/2+1+nx:ny+nx,6)=-1;
    localbd(1+nx+ny:nx/2+nx+ny,7)=1;localbd(nx/2+1+nx+ny:nx+nx+ny,7)=-1;
    localbd(1+2*nx+ny:ny/2+2*nx+ny,8)=1;localbd(ny/2+1+2*nx+ny:ny+2*nx+ny,8)=-1;
    
    localbd(1:nx/4,9)=1;localbd(nx/4+1:nx/2,9)=-1;
    localbd(1+nx:ny/4+nx,10)=1;localbd(ny/4+1+nx:ny/2+nx,10)=-1;
    localbd(1+nx+ny:nx/4+nx+ny,11)=1;localbd(nx/4+1+nx+ny:nx/2+nx+ny,11)=-1;
    localbd(1+nx*2+ny:ny/4+nx*2+ny,12)=1;localbd(nx/4+1+nx*2+ny:ny/2+nx*2+ny,12)=-1;  
    localbd(1+nx/2:3/4*nx,13)=1;localbd(3*nx/4+1:nx,13)=-1;
    localbd(1+ny/2+nx:3*ny/4+nx,14)=1;localbd(3*ny/4+1+nx:ny+nx,14)=-1;
    localbd(1+nx/2+nx+ny:3*nx/4+nx+ny,15)=1;localbd(3*nx/4+1+nx+ny:nx+nx+ny,15)=-1;
    localbd(1+nx/2+nx*2+ny:3*ny/4+nx*2+ny,16)=1;localbd(3*ny/4+1+nx*2+ny:ny+nx*2+ny,16)=-1;  
%     
end
localnx=nx;localny=ny;
            [localbddof]=getboundarydofa(localnx,localny);

   
     localbd1=zeros((localnx+1)*(localny+1),size(localbd,2));
    localbd1(localbddof,:)=localbd;
 