function [ uf,tfineon,Bee,uee]=condense_solve1(A,F,p,bdbasisdof,innerbasisdof,invAii)
% if p==1
%     fineons= tic;
%     uf=A\F;
%     tfineon = toc(fineons);
%     Bee=A;
% else
Aee=A(bdbasisdof,bdbasisdof);Aei=A(bdbasisdof,innerbasisdof);
Aie=A(innerbasisdof,bdbasisdof);
Aii=A(innerbasisdof,innerbasisdof);
Fe=F(bdbasisdof,:);
Fi=F(innerbasisdof,:);
clear A
if nargin==5
[L,U,P,Q] = lu(Aii);
Bee=Aee-Aei*(Q*(U\(L\(P*(Aie)))));
Fee1=Fe-Aei*(Q*(U\(L\(P*(Fi)))));
% Bee=Aee-Aei*(Aii\Aie);
% Fee1=Fe-Aei*(Aii\Fi);
clear L U P Q
% [L,U]=lu(Aii);
% Bee=Aee-Aei*( (U\(L\Aie)));
% Fee1=Fe-Aei*( (U\(L\Fi)));
else %%%% need to preasemble the local inversion matrix
Bee=Aee-Aei*(invAii*Aie);
Fee1=Fe-Aei*(invAii*Fi);
end
clear Aee
fineons= tic;
uee=Bee\Fee1;
tfineon = toc(fineons);

uf=zeros(size(F));
uf(bdbasisdof,:)=uee;
if nargin==5
[L,U,P,Q] = lu(Aii);
uii=Q*(U\(L\(P*(Fi-Aie*uee))));  
% uii=Aii\(Fi-Aie*uee); 
% uii=U\(L\ (Fi-Aie*uee));


else
uii=invAii*(Fi-Aie*uee);
end
uf(innerbasisdof,:)=uii;
% uf=real(uf);

% end