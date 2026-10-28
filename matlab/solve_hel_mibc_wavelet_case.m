function out = solve_hel_mibc_wavelet_case(ctx, ell, overlap_elements, opts)
%SOLVE_HEL_MIBC_WAVELET_CASE Latest-paper edge multiscale coarse solve.
%
% This implementation has no bubble enrichment. Each subdomain contributes
% 4(2^ell+1) mIBC edge-lift basis functions in two dimensions.

if nargin < 4
    opts = struct();
end

verbose = get_opt(opts, 'verbose', 1);
keep_solution = get_opt(opts, 'keep_solution', 0);
keep_basis = get_opt(opts, 'keep_basis', 0);
keep_coarse_matrix = get_opt(opts, 'keep_coarse_matrix', 0);
keep_coarse_factor = get_opt(opts, 'keep_coarse_factor', 0);
pou_type = get_opt(opts, 'pou_type', 'cubic');
dd_info = get_opt(opts, 'dd_info', []);
trace_family = get_opt(opts, 'trace_family', 'independent_edges');

Nx = ctx.Nx;
Ny = ctx.Ny;
nx = ctx.nx;
ny = ctx.ny;
p = ctx.p;
ny_total = ctx.ny_total;
nx_total = ctx.nx_total;
nlocalbasis = 4 * (2^ell + 1);
if strcmp(trace_family,'continuous_boundary'), nlocalbasis = 4*2^ell; end

if overlap_elements < 0 || overlap_elements ~= round(overlap_elements)
    error('overlap_elements must be a nonnegative integer.');
end
if overlap_elements > nx || overlap_elements > ny
    error('The paper assumes delta <= H, hence overlap_elements must not exceed nx or ny.');
end

H = max(ctx.Hx, ctx.Hy);
artificial_ibc_coeff = 1i / H;
ibc_opts = struct('ibc_coeff', artificial_ibc_coeff);

if verbose
    fprintf(['Building mIBC coarse space: ell=%d, basis=%d per patch, ', ...
        'delta/h=%d, POU=%s, artificial IBC=+i/H\n'], ...
        ell, nlocalbasis, overlap_elements, pou_type);
end

total_timer = tic;
basis_timer = tic;
if isempty(dd_info)
    pu1dx = generateoverlappu1d(nx * p, Nx, overlap_elements * p, pou_type);
    pu1dy = generateoverlappu1d(ny * p, Ny, overlap_elements * p, pou_type);
    all_dof = reshape(1:(nx_total * p + 1) * (ny_total * p + 1), ...
        nx_total * p + 1, ny_total * p + 1)';
else
    if dd_info.overelement ~= overlap_elements
        error('dd_info overlap (%d) does not match overlap_elements (%d).', ...
            dd_info.overelement, overlap_elements);
    end
    if ~strcmpi(dd_info.pou_type, pou_type)
        error('dd_info POU type (%s) does not match requested type (%s).', ...
            dd_info.pou_type, pou_type);
    end
end

dofc = Nx * Ny * nlocalbasis;
max_local_nodes = (nx + 2 * overlap_elements + 1) * ...
    (ny + 2 * overlap_elements + 1);
Phi = spalloc(size(ctx.Astiffd, 1), dofc, ...
    Nx * Ny * nlocalbasis * max_local_nodes);
progress_stride = max(1, floor(Nx * Ny / 5));

for iy = 1:Ny
    for ix = 1:Nx
        ie = (iy - 1) * Nx + ix;
        if verbose && mod(Nx * Ny - ie, progress_stride) == 0
            fprintf('remaining coarse patches: %d\n', Nx * Ny - ie);
        end

        localye = max(1, (iy - 1) * ny + 1 - overlap_elements): ...
            min(ny_total, iy * ny + overlap_elements);
        localxe = max(1, (ix - 1) * nx + 1 - overlap_elements): ...
            min(nx_total, ix * nx + overlap_elements);
        if isempty(dd_info)
            localyp = localye(1):(localye(end) + 1);
            localxp = localxe(1):(localxe(end) + 1);
            localny = numel(localye);
            localnx = numel(localxe);
            localdof = all_dof(localyp, localxp);
            localdof = sort(localdof(:));
            global_pu = kron(pu1dy(:, iy), pu1dx(:, ix));
            local_pu = global_pu(localdof);
        else
            localny = dd_info.alllocalny(ie);
            localnx = dd_info.alllocalnx(ie);
            localdof = dd_info.alllocaldof{ie};
            local_pu = dd_info.alllocalpu{ie};
        end

        localF = build_hel_mibc_edge_loads(localnx, localny, ell, ...
            ctx.element_Mass1d, p, trace_family);
        if isempty(dd_info)
            localK = assemble_hel_local_modified_ibc(ctx, localye, localxe, ...
                localnx, localny, ctx.element_Mass1d, ibc_opts);
            basisfun = localK \ localF;
        else
            basisfun = dd_info.alllocalU{ie} \ ...
                (dd_info.alllocalL{ie} \ localF);
        end
        basisfun = basisfun .* local_pu;

        cols = (ie - 1) * nlocalbasis + (1:nlocalbasis);
        Phi(localdof, cols) = sparse(basisfun);
    end
end
basis_time = toc(basis_timer);

Kms = Phi' * ctx.Astiffd * Phi;
Fms = Phi' * ctx.Fd;
coarse_factor_time = 0;
if keep_coarse_factor
    factor_timer = tic;
    [Lc, Uc, Pc, Qc] = lu(Kms);
    coarse_factor_time = toc(factor_timer);
    solve_timer = tic;
    ums = Qc * (Uc \ (Lc \ (Pc * Fms)));
    coarse_solve_time = toc(solve_timer);
else
    solve_timer = tic;
    ums = Kms \ Fms;
    coarse_solve_time = toc(solve_timer);
end
umsfine = Phi * ums;

l2e = NaN; h1e = NaN; hk1e = NaN;
if ~isempty(ctx.u)
    diff = ctx.u - umsfine;
    l2e = relative_matrix_error(diff, ctx.u, ctx.Massweight);
    h1e = relative_matrix_error(diff, ctx.u, ctx.H1matrix);
    hk1e = relative_matrix_error(diff, ctx.u, ctx.Hk1matrix);
end
rel_res = norm(ctx.Astiffd * umsfine - ctx.Fd) / norm(ctx.Fd);

out = struct();
out.ell = ell;
out.trace_family = trace_family;
out.nbasis_per_edge = 2^ell + 1;
out.nbasis_per_patch = nlocalbasis;
out.overlap_elements = overlap_elements;
out.delta = overlap_elements * ctx.hx;
out.pou_type = pou_type;
out.ibc_coeff = artificial_ibc_coeff;
out.dofc = dofc;
out.nfdof = size(ctx.Astiffd, 1);
out.l2e = l2e;
out.h1e = h1e;
out.hk1e = hk1e;
out.rel_res = rel_res;
out.basis_time = basis_time;
out.coarse_factor_time = coarse_factor_time;
out.coarse_solve_time = coarse_solve_time;
out.total_time = toc(total_timer);

if keep_solution
    out.ums = ums;
    out.umsfine = umsfine;
end
if keep_basis
    out.Phi = Phi;
end
if keep_coarse_matrix
    out.Kms = Kms;
    out.Fms = Fms;
    out.phiK = Phi' * ctx.Astiffd;
end
if keep_coarse_factor
    out.Lc = Lc;
    out.Uc = Uc;
    out.Pc = Pc;
    out.Qc = Qc;
end

if verbose
    fprintf(['coarse dof=%d, relative L2/H1/Hk1 errors=%.3e %.3e %.3e, ', ...
        'true coarse residual=%.3e\n'], dofc, l2e, h1e, hk1e, rel_res);
    fprintf('basis/coarse-factor/coarse-solve/total time: %.2f %.2f %.2f %.2f seconds\n', ...
        basis_time, coarse_factor_time, coarse_solve_time, out.total_time);
end
end

function e = relative_matrix_error(diff, ref, W)
e = sqrt(max(real(diff' * W * diff), 0)) / sqrt(max(real(ref' * W * ref), eps));
end

function val = get_opt(opts, name, default_val)
if isfield(opts, name)
    val = opts.(name);
else
    val = default_val;
end
end
