function [localF, meta] = build_hel_mibc_edge_loads(localnx, localny, ell, element_Mass1d, p, trace_family)
%BUILD_HEL_MIBC_EDGE_LOADS Assemble (psi,w)_boundary edge by edge.
%
% Each edge has the hierarchical basis of V_ell: two endpoint functions at
% level zero and the odd-node hats added at every subsequent level. Hence
% there are 2^ell+1 functions per edge and 4(2^ell+1) local basis loads.

if nargin < 5
    p = 1;
end
if nargin < 6, trace_family = 'independent_edges'; end
if strcmp(trace_family, 'continuous_boundary')
    assert(p == 1 && ell >= 0 && ell == round(ell) && ...
        min(localnx,localny)>=2^ell, 'Invalid continuous trace grid or level.');
    ntrace = 4 * 2^ell;
    bd = getboundarydofa(localnx, localny);
    psi = zeros((localnx+1)*(localny+1), ntrace);
    if ell <= 1
        psi(bd,:) = getlocalbdcondition_wavelet_limited(localnx,localny,1,ntrace);
    else
        trace = zeros(numel(bd),ntrace);
        trace(:,1:4) = getlocalbdcondition_wavelet_limited(localnx,localny,1,4);
        lengths = [localnx localny localnx localny];
        offsets = [0 cumsum(lengths(1:3))];
        col=4;
        for level=1:ell
            for edge=1:4
                t=(0:lengths(edge)-1)'/lengths(edge);
                for center=(1:2:2^level-1)/2^level
                    col=col+1;
                    trace(offsets(edge)+(1:lengths(edge)),col)= ...
                        max(1-abs(t-center)*2^level,0);
                end
            end
        end
        psi(bd,:)=trace;
    end
    Mbd = assemble_massedge(ones(localny,localnx),element_Mass1d,p);
    localF = Mbd * psi;
    meta = struct('ell',ell,'nedge_basis',ntrace,'trace_family',trace_family);
    return;
end
assert(strcmp(trace_family,'independent_edges'), 'Unknown trace family.');
if ell < 0 || ell ~= round(ell)
    error('ell must be a nonnegative integer.');
end

nxn = p * localnx + 1;
nyn = p * localny + 1;
nlocaldof = nxn * nyn;
n1d = 2^ell + 1;
nedge_basis = 4 * n1d;

if min(localnx, localny) * p < 2^ell
    error('The local fine-grid trace must have at least 2^ell elements per edge.');
end

unit_coeff = ones(localny, localnx);
[~, Ml, Mr, Mb, Mt] = assemble_massedge(unit_coeff, element_Mass1d, p);
dof = reshape(1:nlocaldof, nxn, nyn)';

Hbottom = hierarchical_trace_basis(linspace(0, 1, nxn)', ell);
Hright = hierarchical_trace_basis(linspace(0, 1, nyn)', ell);
Htop = hierarchical_trace_basis(linspace(0, 1, nxn)', ell);
Hleft = hierarchical_trace_basis(linspace(0, 1, nyn)', ell);

localF = zeros(nlocaldof, nedge_basis);
col = 0;
[localF, col] = add_edge_loads(localF, col, Mb, dof(1, :)', Hbottom);
[localF, col] = add_edge_loads(localF, col, Mr, dof(:, end), Hright);
[localF, col] = add_edge_loads(localF, col, Mt, dof(end, :)', Htop);
[localF, ~] = add_edge_loads(localF, col, Ml, dof(:, 1), Hleft);

meta = struct();
meta.ell = ell;
meta.n1d = n1d;
meta.nedge_basis = nedge_basis;
meta.edge_order = {'bottom', 'right', 'top', 'left'};
end

function [localF, col] = add_edge_loads(localF, col, Medge, edgedof, H)
nlocaldof = size(localF, 1);
nedge = size(H, 2);
psi = zeros(nlocaldof, nedge);
psi(edgedof, :) = H;
cols = col + (1:nedge);
localF(:, cols) = Medge * psi;
col = col + nedge;
end

function H = hierarchical_trace_basis(t, ell)
% Ordering: endpoint functions first, then odd dyadic hats level by level.
H = [1 - t, t];
for level = 1:ell
    centers = (1:2:(2^level - 1)) / 2^level;
    for center = centers
        H(:, end + 1) = max(1 - abs(t - center) * 2^level, 0); %#ok<AGROW>
    end
end
end
