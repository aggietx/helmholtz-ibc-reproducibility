[allbddof_condense,allinnerdof_condense]=sortdof1(Nx,Ny,nx,ny,p);
Minvii=computeinnerinver_2d_per(Nx*Ny,K(allinnerdof_condense,allinnerdof_condense),...
    p,perindex,local_permodel);
u=condense_solve1(K,F,p,allbddof_condense,allinnerdof_condense,Minvii);
u=real(u);
return
Ac=sparse(ixc,iyc,localmsmatrix);
[allbddofc,allinnerdofc]=sortdof(Nx,Ny,p);
invAcinner=computeinnerinver(Nx*Ny,Ac(allinnerdofc,allinnerdofc),p,2);
tmsoff=toc(msoffs);
% clear ixc iyc localmsstiffmatrix
%% solve coarse
disp('solve coarse system...')
% [Ac,Fc]=diri_set(Ac,Fc,bddofc,0);
% % Ac=Ac+10^(-8)*diag(diag(Ac));

[uc,tmson,Acee,uee]=condense_solve1(Ac,Fc,p,allbddofc,allinnerdofc,invAcinner);
% % msons=tic;uee=Ac(allbddofc,allbddofc)\Fc(allbddofc,:);
% % uc=zeros(size(Fc));uc(allbddofc,:)=uee;tmson=toc(msons);
% msons=tic;uc=real(Ac\Fc);tmson=toc(msons);
% figure();imagesc(reshape(uc,p*Nx+1,p*Ny+1)');colorbar;title('MsFEM  solution (coarse)');
% return

%% recover full fine scale solution
% umsfine=recoverfinep(Ny,Nx,nx0,ny0,trialbasisfull,uc,node_matrixc,p);
umsfine=recover_meshp_elliptic(Ny,Nx,nx0,ny0,trailbasis_meshnodes,uc,node_matrixc,p);

umsfine=real(umsfine);
umsm=reshape(umsfine,1*nx_total+1,1*ny_total+1)';
umsfine0 =umsm(damp_bottom*1+1:end-damp_top*1,damp_left*1+1:end-damp_right*1);

if solve_fine==1
figure();imagesc([0,Ly_true],[0,Lx_true],ufine0);colorbar;title('Sem solution');
end
% if samplingratio==6
% load u6a;uref=u6a; %%%% top source
% % load u6;uref=u6; %% ceneter source
% elseif samplingratio==10
% load u10a;uref=u10a; %%%% top source
% % load u10;uref=u10;%% ceneter source
% end 
load uref
% load ufine0;uref=ufine0;
e=uref-umsfine0;ea=e.';ua=uref.';
l2error=sqrt(ea(:).'*ea(:))/sqrt(ua(:).'*ua(:));
fprintf('Relative L2 ERROR is  %2.4e\n',l2error);

return
neighorbood_dof=reshape(1:(2*nx+1)*(2*ny+1),2*nx+1,2*ny+1)';

left_bottomdof=neighorbood_dof(1:ny+1, 1:nx+1);left_bottomdof=sort(left_bottomdof(:));
right_bottomdof=neighorbood_dof(1:ny+1, nx+1:end);right_bottomdof=sort(right_bottomdof(:));
left_topdof=neighorbood_dof(ny+1:end,  1:nx+1);left_topdof=sort(left_topdof(:));
right_topdof=neighorbood_dof(ny+1:end, nx+1:end);right_topdof=sort(right_topdof);
iivalue=ones((2*nx+1)*(2*ny+1),1);
iy=2;ix=2;
 
                  ie=(iy-1)*(Nx+1)+ix;
      if mod((Nx+1)*(Ny+1)-ie,floor((Nx+1)*(Ny+1)/5))==0
          fprintf('rest coarse neighborhood is %d\n',(Nx+1)*(Ny+1)-ie);
      end
      
for kk=1:4
    if kk==1 && elt(kk)>0
     %---bottom left element.
     jj=ceil(elt(kk)/Nx);
     ii=elt(kk)-(jj-1)*Nx;
     gi_element=getCoarseElemIndex(ii, jj, Nx, Ny, nx, ny);%%%% dof in a elements
%      phimatrix(nind, gi_element) = all_nodalbasis(:,4,  elt(kk)); 
%     iix(left_bottomdof,ie)=gi_element;
%     iivalue(left_bottomdof)=polybasis(:,4); 
    iivalue(left_bottomdof)=alltrialbasis(:,4,  elt(kk)); 
    elseif kk==2 && elt(kk)>0
    %---bottom right element.  
      jj=ceil(elt(kk)/Nx);
      ii=elt(kk)-(jj-1)*Nx;
      gi_element=getCoarseElemIndex(ii, jj, Nx, Ny, nx, ny);
%       phimatrix(nind, gi_element) = all_nodalbasis(:,3, elt(kk));
      
%     iix(right_bottomdof,ie)=gi_element;
%     iivalue(right_bottomdof)=polybasis(:,3);  
    iivalue(right_bottomdof)=alltrialbasis(:,3,  elt(kk)); 
    elseif kk==3 && elt(kk)>0
    %---top right element.   
      jj=ceil(elt(kk)/Nx);
      ii=elt(kk)-(jj-1)*Nx;
      gi_element=getCoarseElemIndex(ii, jj, Nx, Ny, nx, ny);
%      phimatrix(nind, gi_element) = all_nodalbasis(:,1,  elt(kk));
     
%          iix(right_topdof,ie)=gi_element;
%     iivalue(right_topdof)=polybasis(:,1);
    iivalue(right_topdof)=alltrialbasis(:,1,  elt(kk)); 
    elseif kk==4 && elt(kk)>0 
    %---top left element. 
      jj=ceil(elt(kk)/Nx);
      ii=elt(kk)-(jj-1)*Nx;
      gi_element=getCoarseElemIndex(ii, jj, Nx, Ny, nx, ny);
%      phimatrix(nind, gi_element) = all_nodalbasis(:,2,  elt(kk));
    
%     iix(left_topdof,ie)=gi_element;
%     iivalue(left_topdof)=polybasis(:,2);  
     iivalue(left_topdof)=alltrialbasis(:,2,  elt(kk)); 
    end
  end
% reshape(iivalue,2*nx+1,2*ny+1)';
imagescsquare(iivalue)
return


i=10000;
    M1 = rand (i,i);
    M2 = rand (i,i); 
    tic
    M = M1 * M2;
    toc
    M1 = distributed(M1);
    M2 = distributed(M2);
    tic
    M = M1 * M2;
    toc
    return
close all

[xo,yo]=getcoor_sem(-1,1,-1,1,nx0,ny0,p);
[xn,yn]=getcoor_sem_a(-1,1,-1,1,nx0,ny0,p);
polybasiso=get_coarseleg1(xo,yo,p);
polybasisn=get_coarseleg1(xn,yn,p);
norm(polybasiso-polybasisn)/norm(polybasiso)
return
for i=1:(p+1)^2
    meshsquare(polybasiso(:,i));title('old')
    pause()
    meshsquare(polybasisn(:,i));title('new')
    pause()
    close all
end
    
return
fid=fopen('frac_basis3.bin','w+');fwrite(fid,frac_basis3,'single');fclose(fid); %%%% save
return
[Astiff,F]=diri_set(Astiff,F,bdnodes,0);
u=Astiff\F;
return

p=3;
[x,y]=getcoor_sem(-1,1,-1,1,10,10,p);
ndof=size(x,1);
a = get_quadrature_point_weight(p+1,-1,1); %% a is quad points 
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

return
close all
i=5;
% load seif
% i=4;load seifq4
% load seicq4
refsei=seif4;
seir=seif;
plot(refsei(i,:),'r');hold on;plot(seir(i,:),'g');%plot(seicq4(i,:),'b');
norm(seir-refsei)/norm(refsei)
return
close all;
for i=1:size(a,2)
     meshsquare(a(:,i));title('poly'); pause();close all
      meshsquare(tempmsbasis(:,i));title('ms');pause();close all
end
for i=1:36
    norm(msbasis(:,:,i)-a(:,:,i))
end

for i=1:9
%        meshsquare(a(:,i,1));title('poly'); pause();close all
%       meshsquare(msbasis(:,i,1));title('ms');pause();close all  
      imagescsquare(a(:,i,1));title('poly'); pause();close all
     imagescsquare(msbasis(:,i,1));title('ms');pause();close all  
      norm(a(:,i,1)-msbasis(:,i,1))
end
return
n_node=7;
[x,wx]=xwlgl(n_node,0,1);
[y,wy]=xwlgl(n_node,0,1);
deri_x=derlgl(x,n_node);%first derivative GLL matrix % % saved in each column 
deri_y=derlgl(y,n_node);
jacx=hx/2;
jacy=hy/2;
% % jac=jacx*jacy;

npdx=length(wx);
npdy=length(wy);
jacyx=jacy/jacx;
jacxy=jacx/jacy;
mn=npdx*npdy;
[wx1,wy1]=meshgrid(wx,wy); ww=wx1.*wy1; ww=ww'; ww=ww(:); clear wx1 wy1;
coeff=1*rand(size(ww,1),1);
sum(ww)
% ww'*coeff/4


return


localA=assemble_weight_stiff_sem(local_Astiff,local_stiff_coeff_ms,localix,localiy);localA0=localA;
    localMass=assemblemassq1(local_Mass,local_mass_coeff_ms,localix,localiy);
    Hx=nx0*hx;Hy=ny0*hy;nlocaldof=(2*nx0+1)*(2*ny0+1);
%     F=zeros((2*nx0+1)*(2*ny0+1),9);
%% step 1(boundary nodes)
    x=x1/Hx;y=y1/Hy;x=x.';y=y.';x=x(:);y=y(:);
    F1=(1-2*x).*(1-x).*(1-2*y).*(1-y);
    F2=4.*x.*(1-x).*(1-2*y).*(1-y);
    F3=x.*(2*x-1).*(1-2*y).*(1-y);
    
    F4=(1-2*x).*(1-x).*4.*y.*(1-y);
    F5=4.*x.*(1-x).*4.*y.*(1-y);
    F6=x.*(2*x-1).*4.*y.*(1-y);
    
    F7=(1-2*x).*(1-x).*y.*(2*y-1);
    F8=4.*x.*(1-x).*y.*(2*y-1);
    F9=x.*(2*x-1).*y.*(2*y-1);
    F=[F1,F2,F3,F4,F5,F6,F7,F8,F9];
    F(local_freedofa,:)=0;
  localA(local_bddofa,:)=sparse(1:size(local_bddofa,1),local_bddofa',...
        ones(1,size(local_bddofa,1)),size(local_bddofa,1),(2*nx0+1)*(2*ny0+1));
msbasistemp= localA\F;
%% step 2 (interior nodes)
forceinner=32*( x.*(1-x)+y.*(1-y));
Fi=localMass*forceinner;Fi(local_bddofa,:)=0;
% Fi=ones((2*nx0+1)*(2*ny0+1),1); Fi(local_bddofa,:)=0;
basisi=localA\Fi;
 meshsquare(basisi)
 meshsquare(4.*x.*(1-x).*4.*y.*(1-y))
return
%% combine

invalue=basisi( (size(basisi,1)-1)/2+1);
basisirep=repmat(basisi,1,9);
coef=-msbasistemp((size(basisi,1)-1)/2+1,:)/invalue;
coefrep=repmat(coef,nlocaldof,1);
msbasis=msbasistemp+coefrep.*basisirep;
msbasis(:,5)=basisi/invalue;
