function Minvii=computeinnerinver(Nblock,Mii,p,d)
%%% Nblock=Nx*Ny;
%%% d is the dimension of the problem
ninndof=(p-1)^d;
localMinvii=cell(Nblock,1);
local_Mii=cell(Nblock,1);

for ii=1:Nblock
    localdof= (ii-1)*ninndof+1:ii*ninndof;
    local_Mii{ii}=Mii(localdof,localdof);
end
parfor ii=1:Nblock
%     localdof= (ii-1)*ninndof+1:ii*ninndof;
%     localMii=Mii(localdof,localdof);
    localMii=local_Mii{ii};
    localMinvii{ii}=localMii\speye(ninndof);
end
Minvii=blkdiag(localMinvii{:});
