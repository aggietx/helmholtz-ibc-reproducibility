function [quadpoint,quadweight] = gllquad(p,a,b) 
%%%% the gll points and weight in the interval [a,b],default interval is
%%%% [-1,1]
%%%% example:
%%%% [x,weight] = get_gll_weight_position(4);
np=p+1;
if np <= 1 
  es='The number of the quadrature nodes should be greater than 1';
  disp(es); quadpoint = []; quadweight = [];
  return 
  
elseif np==2
    quadpoint=[-1;1];quadweight=[1;1];
    
elseif np==3
     quadpoint=[-1;0;1];quadweight=[0.333333333333333;
   1.333333333333333;
   0.333333333333333];
elseif np==4
  quadpoint=[-1.000000000000000;
  -0.447213595499958;
   0.447213595499958;
   1.000000000000000]; 
quadweight= [0.166666666666667;
   0.833333333333334;
   0.833333333333334;
   0.166666666666667];
elseif np==5
  quadpoint=[-1.000000000000000;
  -0.654653670707977;
   0.000000000000000;
   0.654653670707977;
   1.000000000000000];
quadweight=[ 0.100000000000000;
   0.544444444444444;
   0.711111111111111;
   0.544444444444444;
   0.100000000000000];
elseif np==6
 quadpoint=[-1.000000000000000;
  -0.765055323929465;
  -0.285231516480645;
   0.285231516480645;
   0.765055323929465;
   1.000000000000000];
quadweight=[0.066666666666667;
   0.378474956297847;
   0.554858377035486;
   0.554858377035486;
   0.378474956297847;
   0.066666666666667];
elseif np==7
 quadpoint=[-1.000000000000000;
  -0.830223896278567;
  -0.468848793470714;
                   0;
   0.468848793470714;
   0.830223896278567;
   1.000000000000000];  
quadweight=[ 0.047619047619048;
   0.276826047361566;
   0.431745381209863;
   0.487619047619048;
   0.431745381209863;
   0.276826047361566;
   0.047619047619048];
else
quadpoint=zeros(np,1);
% w=zeros(np,1);
n=np-1;
quadpoint(1)=-1;quadpoint(np)=1;
x1=jacobi_roots(n-1,1,1);
quadpoint(2:n)=x1;
coef=2/(n*np);
quadweight=coef./(pnleg(quadpoint,n).^2);
end

if nargin ==3 
   bma=(b-a)*.5;
   bpa=(b+a)*.5;
   quadpoint=bma*quadpoint+bpa;
   quadweight=quadweight*bma;
end

return
