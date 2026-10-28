function Minvii=computeinnerinver_2d_per(Nblock,Mii,p,perindex,local_permodel)
%%% actually not only for q1
[ny,nx]=size(local_permodel);
% localnonzerodof=getnonzerodof(nx,ny,local_permodel,getnod(p,nx,ny));
per_index=find(perindex==1);

ninndof=(nx*p-1)*(ny*p-1);
localMinvii=cell(Nblock,1);
local_Mii=cell(Nblock,1);

for ii=1:Nblock
    localdof= (ii-1)*ninndof+1:ii*ninndof;
    local_Mii{ii}=Mii(localdof,localdof);
end
sampleMii=local_Mii{per_index(1)};
localnonzerodof=[];
for i=1:size(sampleMii,1)
    if sampleMii(i,i)~=0
    localnonzerodof=[localnonzerodof,i];
    end
end
% whos localnonzerodof
parfor ii=1:Nblock
% for ii=1:Nblock    
    if ismember(ii,per_index)
        localMii=local_Mii{ii};
%         localM=sparse(1:ninndof,1:ninndof,zeros(ninndof,1));
        localM=zeros(ninndof,ninndof);
        localM(localnonzerodof,localnonzerodof)=localMii(localnonzerodof,localnonzerodof)\speye(length(localnonzerodof));
         localMinvii{ii}=localM;
    else
        %     localdof= (ii-1)*ninndof+1:ii*ninndof;
%     localMii=Mii(localdof,localdof);
    localMii=local_Mii{ii};
    localMinvii{ii}=localMii\speye(ninndof);

    end
end
Minvii=blkdiag(localMinvii{:});
