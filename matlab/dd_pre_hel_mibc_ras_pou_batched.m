function z = dd_pre_hel_mibc_ras_pou_batched(dd, inc)
%DD_PRE_HEL_MIBC_RAS_POU_BATCHED Apply weighted RAS using batched local solves.

if ~isfield(dd, 'factor_groups') || isempty(dd.factor_groups)
    z = dd_pre_hel_mibc_ras_pou(dd.alllocalL, dd.alllocalU, ...
        dd.alllocaldof, dd.alllocalpu, inc);
    return;
end

batch_size = dd.batch_size;
z = zeros(size(inc), 'like', inc);
for ig = 1:numel(dd.factor_groups)
    patches = dd.factor_groups{ig};
    if isempty(patches)
        continue;
    end
    L = dd.alllocalL{patches(1)};
    U = dd.alllocalU{patches(1)};
    nlocal = numel(dd.alllocaldof{patches(1)});

    for first = 1:batch_size:numel(patches)
        ids = patches(first:min(first + batch_size - 1, numel(patches)));
        rhs = zeros(nlocal, numel(ids), 'like', inc);
        for j = 1:numel(ids)
            rhs(:, j) = inc(dd.alllocaldof{ids(j)});
        end
        localz = U \ (L \ rhs);
        for j = 1:numel(ids)
            ie = ids(j);
            dof = dd.alllocaldof{ie};
            z(dof) = z(dof) + dd.alllocalpu{ie} .* localz(:, j);
        end
    end
end
end
