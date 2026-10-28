function z = dd_pre_hel_mibc_ras_pou(alllocalL, alllocalU, ...
    alllocaldof, alllocalpu, inc)
%DD_PRE_HEL_MIBC_RAS_POU Apply the latest-paper weighted RAS-imp operator.

z = zeros(size(inc));
for ie = 1:numel(alllocaldof)
    localdof = alllocaldof{ie};
    localz = alllocalU{ie} \ (alllocalL{ie} \ inc(localdof));
    z(localdof) = z(localdof) + alllocalpu{ie} .* localz;
end
end
