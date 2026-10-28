function info = build_hel_dd_info_modified_ibc(ctx, overelement, opts)
%BUILD_HEL_DD_INFO_MODIFIED_IBC Prepare RAS local LU factors for Helmholtz.
%
% The local DD problems use physical ABC on the global boundary and the
% modified impedance coefficient on artificial subdomain boundaries.

if nargin < 3
    opts = struct();
end

verbose = get_opt(opts, 'verbose', 1);
tao = get_opt(opts, 'tao', 0);
pou_type = get_opt(opts, 'pou_type', 'cubic');
reuse_homogeneous_factors = get_opt(opts, 'reuse_homogeneous_factors', false);

Nx = ctx.Nx;
Ny = ctx.Ny;
nx = ctx.nx;
ny = ctx.ny;
p = ctx.p;
ny_total = Ny * ny;
nx_total = Nx * nx;

if isfield(ctx, 'element_Mass1d')
    element_Mass1d = ctx.element_Mass1d;
else
    [element_Mass1d, ~] = get_elementmatrix_blend1d(ctx.hx, tao, p);
end

all_dof = reshape(1:(nx_total * p + 1) * (ny_total * p + 1), ...
    nx_total * p + 1, ny_total * p + 1).';

alllocaldof = cell(Nx * Ny, 1);
alllocaldof0 = cell(Nx * Ny, 1);
all_truedof = cell(Nx * Ny, 1);
alllocalL = cell(Nx * Ny, 1);
alllocalU = cell(Nx * Ny, 1);
allfactor_group = zeros(Nx * Ny, 1);
alllocalny = zeros(Nx * Ny, 1);
alllocalnx = zeros(Nx * Ny, 1);
localye_cell = cell(Nx * Ny, 1);
localxe_cell = cell(Nx * Ny, 1);
alllocalpu = cell(Nx * Ny, 1);
artificial_ibc_coeff = NaN;

if reuse_homogeneous_factors && ~is_constant_context(ctx)
    warning('Factor reuse was disabled because the coefficient fields are not constant.');
    reuse_homogeneous_factors = false;
end

if verbose
    fprintf('preparing modified-IBC RAS information, dd overlap delta/h=%d\n', overelement);
end

setup_timer = tic;
for iy = 1:Ny
    for ix = 1:Nx
        ie = (iy - 1) * Nx + ix;

        localye = max(1, (iy - 1) * ny + 1 - overelement): ...
            min(ny_total, iy * ny + overelement);
        localxe = max(1, (ix - 1) * nx + 1 - overelement): ...
            min(nx_total, ix * nx + overelement);
        localyp = max(1, (iy - 1) * ny * p + 1 - overelement * p): ...
            min(ny_total * p + 1, iy * ny * p + overelement * p + 1);
        localxp = max(1, (ix - 1) * nx * p + 1 - overelement * p): ...
            min(nx_total * p + 1, ix * nx * p + overelement * p + 1);

        localdof = all_dof(localyp, localxp);
        localdof = sort(localdof(:));
        alllocaldof{ie} = localdof;
        alllocalny(ie) = numel(localye);
        alllocalnx(ie) = numel(localxe);
        localye_cell{ie} = localye;
        localxe_cell{ie} = localxe;

        localyp0 = max(1, (iy - 1) * ny * p + 1): ...
            min(ny_total * p + 1, iy * ny * p + 1);
        localxp0 = max(1, (ix - 1) * nx * p + 1): ...
            min(nx_total * p + 1, ix * nx * p + 1);
        localdof0 = all_dof(localyp0, localxp0);
        localdof0 = sort(localdof0(:));
        alllocaldof0{ie} = localdof0;

        [tf, true_dof] = ismember(localdof0, localdof);
        if any(~tf)
            error('Core dof not contained in overlapping local dof for block (%d,%d).', ix, iy);
        end
        all_truedof{ie} = true_dof;
    end
end

pu1dx = generateoverlappu1d(nx * p, Nx, overelement * p, pou_type);
pu1dy = generateoverlappu1d(ny * p, Ny, overelement * p, pou_type);
pou_sum = zeros((nx_total * p + 1) * (ny_total * p + 1), 1);
for iy = 1:Ny
    for ix = 1:Nx
        ie = (iy - 1) * Nx + ix;
        global_pu = kron(pu1dy(:, iy), pu1dx(:, ix));
        localdof = alllocaldof{ie};
        alllocalpu{ie} = global_pu(localdof);
        pou_sum(localdof) = pou_sum(localdof) + alllocalpu{ie};
    end
end
pou_defect = norm(pou_sum - 1, inf);
if pou_defect > 1e-12
    error('Full-overlap DD partition of unity defect is %.3e.', pou_defect);
end

progress_stride = max(1, floor(Nx * Ny / 5));
factor_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
factor_count = 0;
for ie = 1:Nx * Ny
    if verbose && mod(Nx * Ny - ie, progress_stride) == 0
        fprintf('rest local DD factors is %d\n', Nx * Ny - ie);
    end

    localny = alllocalny(ie);
    localnx = alllocalnx(ie);
    localye = localye_cell{ie};
    localxe = localxe_cell{ie};
    cache_key = sprintf('%d_%d_%d_%d_%d_%d', localnx, localny, ...
        localxe(1) == 1, localxe(end) == nx_total, ...
        localye(1) == 1, localye(end) == ny_total);
    if reuse_homogeneous_factors && isKey(factor_cache, cache_key)
        factors = factor_cache(cache_key);
        L = factors.L;
        U = factors.U;
        factor_group = factors.group;
    else
        [localK, ~, artificial_ibc_coeff] = assemble_hel_local_modified_ibc(ctx, ...
            localye, localxe, localnx, localny, element_Mass1d, opts);
        [L, U] = lu(localK);
        factor_count = factor_count + 1;
        factor_group = factor_count;
        if reuse_homogeneous_factors
            factor_cache(cache_key) = struct('L', L, 'U', U, ...
                'group', factor_group);
        end
    end
    alllocalL{ie} = L;
    alllocalU{ie} = U;
    allfactor_group(ie) = factor_group;
end
setup_time = toc(setup_timer);

info = struct();
info.overelement = overelement;
info.alllocalL = alllocalL;
info.alllocalU = alllocalU;
info.allfactor_group = allfactor_group;
info.alllocaldof = alllocaldof;
info.all_truedof = all_truedof;
info.alllocaldof0 = alllocaldof0;
info.alllocalpu = alllocalpu;
info.alllocalnx = alllocalnx;
info.alllocalny = alllocalny;
info.localxe_cell = localxe_cell;
info.localye_cell = localye_cell;
info.artificial_ibc_coeff = artificial_ibc_coeff;
info.setup_time = setup_time;
info.local_bc_type = 'modified_ibc';
info.pou_type = pou_type;
info.pou_defect = pou_defect;
info.reuse_homogeneous_factors = reuse_homogeneous_factors;
info.unique_factor_count = factor_count;
info.factor_groups = cell(factor_count, 1);
for ig = 1:factor_count
    info.factor_groups{ig} = find(allfactor_group == ig);
end
info.batch_size = 64;

if verbose
    fprintf(['modified-IBC RAS setup time is %.2f seconds, artificial IBC=%g%+gi, ', ...
        'full-overlap %s POU defect=%.3e, unique local factors=%d\n'], setup_time, ...
        real(artificial_ibc_coeff), imag(artificial_ibc_coeff), pou_type, pou_defect, ...
        factor_count);
end
end

function tf = is_constant_context(ctx)
tf = field_is_constant(ctx.stiff_coeff) && field_is_constant(ctx.mass_coeff) && ...
    field_is_constant(ctx.abc_coeff);
end

function tf = field_is_constant(values)
values = values(:);
scale = max(1, max(abs(values)));
tf = max(abs(values - values(1))) <= 100 * eps(scale);
end

function val = get_opt(opts, name, default_val)
if isfield(opts, name)
    val = opts.(name);
else
    val = default_val;
end
end
