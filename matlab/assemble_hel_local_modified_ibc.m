function [localK, local_edge_mass, artificial_ibc_coeff] = assemble_hel_local_modified_ibc(ctx, ...
    localye, localxe, localnx, localny, element_Mass1d, opts)
%ASSEMBLE_HEL_LOCAL_MODIFIED_IBC Assemble a Helmholtz local operator.
%
% Physical outer boundaries inherit ctx.abc_coeff. Artificial subdomain
% boundaries use the modified impedance coefficient with scale H^{-1}.

if nargin < 7
    opts = struct();
end

p = ctx.p;
if isfield(opts, 'ibc_coeff')
    artificial_ibc_coeff = opts.ibc_coeff;
else
    Hlocal = get_opt(opts, 'ibc_length_scale', max(ctx.Hx, ctx.Hy));
    artificial_wavenumber = get_opt(opts, 'ibc_wavenumber', 1 / Hlocal);
    if isfield(ctx, 'mass_coeff') && isfield(ctx, 'abc_coeff')
        sign_coeff = mean(ctx.abc_coeff(:) ./ sqrt(ctx.mass_coeff(:)));
    else
        sign_coeff = -1i;
    end
    artificial_ibc_coeff = sign_coeff * artificial_wavenumber;
end

[iix0, iiy0] = get_ix_iy0(p, localnx, localny);
localA = assemblymatrix(iix0, iiy0, ctx.stiff_coeff(localye, localxe), ctx.element_Astiff);
localM = assemblymatrix(iix0, iiy0, ctx.mass_coeff(localye, localxe), ctx.element_Mass);

coeff_left = zeros(localny, localnx);
coeff_right = zeros(localny, localnx);
coeff_bottom = zeros(localny, localnx);
coeff_top = zeros(localny, localnx);

nx_total = ctx.nx_total;
ny_total = ctx.ny_total;

if localxe(1) == 1
    coeff_left(:, 1) = ctx.abc_coeff(localye, 1);
else
    coeff_left(:, 1) = artificial_ibc_coeff;
end

if localxe(end) == nx_total
    coeff_right(:, end) = ctx.abc_coeff(localye, nx_total);
else
    coeff_right(:, end) = artificial_ibc_coeff;
end

if localye(1) == 1
    coeff_bottom(1, :) = ctx.abc_coeff(1, localxe);
else
    coeff_bottom(1, :) = artificial_ibc_coeff;
end

if localye(end) == ny_total
    coeff_top(end, :) = ctx.abc_coeff(ny_total, localxe);
else
    coeff_top(end, :) = artificial_ibc_coeff;
end

[~, Ml, ~, ~, ~] = assemble_massedge(coeff_left, element_Mass1d, p);
[~, ~, Mr, ~, ~] = assemble_massedge(coeff_right, element_Mass1d, p);
[~, ~, ~, Mb, ~] = assemble_massedge(coeff_bottom, element_Mass1d, p);
[~, ~, ~, ~, Mt] = assemble_massedge(coeff_top, element_Mass1d, p);

localK = localA - localM + Ml + Mr + Mb + Mt;

if nargout >= 2
    unit_coeff = ones(localny, localnx);
    [local_edge_mass, ~, ~, ~, ~] = assemble_massedge(unit_coeff, element_Mass1d, p);
end
end

function val = get_opt(opts, name, default_val)
if isfield(opts, name)
    val = opts.(name);
else
    val = default_val;
end
end
