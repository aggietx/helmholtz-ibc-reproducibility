function ctx = build_hel_mibc_fine_context(opts)
%BUILD_HEL_MIBC_FINE_CONTEXT Assemble a parameterized 2-D Helmholtz problem.
%
% The sign convention follows the latest mIBC paper:
%   -Delta u - k^2 u = f,  partial_n u + i k u = 0.

% The fine and edge matrices use the same blended Q1 quadrature.

if nargin < 1
    opts = struct();
end

L = get_opt(opts, 'domain_length', 1);
fine_n = get_opt(opts, 'fine_n', 32);
coarse_n = get_opt(opts, 'coarse_n', 4);
k = get_opt(opts, 'k', 4);
p = get_opt(opts, 'p', 1);
tao = get_opt(opts, 'tao', 0.5);
source_x = get_opt(opts, 'source_x', 0.5 * L);
source_y = get_opt(opts, 'source_y', 0.5 * L);
source_width_cells = get_opt(opts, 'source_width_cells', 1);
verbose = get_opt(opts, 'verbose', 1);

if p ~= 1
    error('The latest-paper 2-D test currently supports Q1 elements only.');
end
if fine_n <= 0 || coarse_n <= 0 || fine_n ~= round(fine_n) || coarse_n ~= round(coarse_n)
    error('fine_n and coarse_n must be positive integers.');
end
if mod(fine_n, coarse_n) ~= 0
    error('fine_n=%d must be divisible by coarse_n=%d.', fine_n, coarse_n);
end
if k <= 0 || L <= 0
    error('k and domain_length must be positive.');
end

Nx = coarse_n;
Ny = coarse_n;
nx = fine_n / coarse_n;
ny = nx;
nx_total = fine_n;
ny_total = fine_n;
hx = L / fine_n;
hy = hx;
Hx = L / coarse_n;
Hy = Hx;

[~, x, y, nod] = mesh_sem_new(L, L, nx_total, ny_total, 1);
[element_Astiff, ~, ~, element_Mass, element_Mass1d] = ...
    get_element_blend_matrix(hx, hy, p, tao);
[ix, iy] = get_ix_iy(p, nod');

stiff_coeff = ones(ny_total, nx_total);
mass_coeff = k^2 * ones(ny_total, nx_total);
abc_coeff = 1i * k * ones(ny_total, nx_total);

Astiff = assemblymatrix(ix, iy, stiff_coeff, element_Astiff);
Massweight = assemblymatrix(ix, iy, ones(ny_total, nx_total), element_Mass);
M = assemblymatrix(ix, iy, mass_coeff, element_Mass);
M_boundary = assemble_massedge(abc_coeff, element_Mass1d, p);
A = Astiff - M + M_boundary;

source_support = source_width_cells * hx;
force = exp(-((x - source_x).^2 + (y - source_y).^2) / source_support^2) / source_support^2;
F = Massweight * force;

if verbose
    fprintf(['Assembled latest-paper 2-D problem: fine=%d x %d, coarse=%d x %d, ', ...
        'k=%g, h=%g, H=%g, ppw=%g, blended tao=%g\n'], ...
        fine_n, fine_n, coarse_n, coarse_n, k, hx, Hx, 2 * pi / (k * hx), tao);
end

fine_timer = tic;
if get_opt(opts, 'solve_reference', true)
    u = A \ F;
else
    u = [];
end
fine_solve_time = toc(fine_timer);
true_rel_res = NaN;
if ~isempty(u), true_rel_res = norm(A * u - F) / norm(F); end

ctx = struct();
ctx.Nx = Nx;
ctx.Ny = Ny;
ctx.nx = nx;
ctx.ny = ny;
ctx.p = p;
ctx.Astiffd = A;
ctx.Fd = F;
ctx.u = u;
ctx.Mass = Massweight;
ctx.Massweight = Massweight;
ctx.H1matrix = Astiff + Massweight;
ctx.H1seminorm_matrix = Astiff;
ctx.Hk1matrix = Astiff + M;
ctx.element_Astiff = element_Astiff;
ctx.element_Mass = element_Mass;
ctx.element_Mass1d = element_Mass1d;
ctx.ixfine = ix;
ctx.iyfine = iy;
ctx.stiff_coeff = stiff_coeff;
ctx.mass_coeff = mass_coeff;
ctx.abc_coeff = abc_coeff;
ctx.hx = hx;
ctx.hy = hy;
ctx.Hx = Hx;
ctx.Hy = Hy;
ctx.nx_total = nx_total;
ctx.ny_total = ny_total;
ctx.nfdof = size(A, 1);
ctx.k0 = k;
ctx.k = k;
ctx.ppw = 2 * pi / (k * hx);
ctx.tao = tao;
ctx.domain_length = L;
ctx.source_x = source_x;
ctx.source_y = source_y;
ctx.source_support = source_support;
ctx.fine_solve_time = fine_solve_time;
ctx.fine_true_rel_res = true_rel_res;

if verbose
    fprintf('Fine direct solve: %.2f seconds, true |Au-f|/|f|=%.3e\n', ...
        fine_solve_time, true_rel_res);
end
end

function val = get_opt(opts, name, default_val)
if isfield(opts, name)
    val = opts.(name);
else
    val = default_val;
end
end
