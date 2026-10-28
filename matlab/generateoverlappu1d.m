function pu1dx = generateoverlappu1d(nx, Nx, overe, pou_type)
%GENERATEOVERLAPPU1D Build 1-D overlap partition-of-unity weights.
%   pou_type choices:
%     'legacy'  - original generatelinearpu1d; requires 2*overe <= nx.
%     'linear'  - normalized linear ramp.
%     'cubic'   - normalized C1 smoothstep ramp.
%     'quintic' - normalized C2 smoothstep ramp.
%     'cosine'  - normalized cosine ramp.

if nargin < 4 || isempty(pou_type)
    pou_type = 'linear';
end
if nargin < 3
    error('generateoverlappu1d requires nx, Nx, and overe.');
end
if nx <= 0 || Nx <= 0
    error('nx and Nx must be positive.');
end
if overe < 0 || abs(overe - round(overe)) > 0
    error('overe must be a nonnegative integer.');
end

pou_type = lower(strtrim(pou_type));
if strcmp(pou_type, 'normalized-linear')
    pou_type = 'linear';
end
if strcmp(pou_type, 'smooth')
    pou_type = 'cubic';
end

if strcmp(pou_type, 'legacy')
    if 2 * overe > nx
        error('legacy POU requires 2*overe <= nx. Use linear/cubic/quintic/cosine for wider overlaps.');
    end
    pu1dx = generatelinearpu1d(nx, Nx, overe);
    return;
end

nnx = nx * Nx + 1;
nodes = (1:nnx)';
pu1dx = zeros(nnx, Nx);

for i = 1:Nx
    left_node = (i - 1) * nx + 1;
    right_node = i * nx + 1;

    if overe == 0
        weight = double(nodes >= left_node & nodes <= right_node);
    else
        left_t = (nodes - (left_node - overe)) / (2 * overe);
        right_t = (right_node + overe - nodes) / (2 * overe);
        left_weight = overlap_transition(left_t, pou_type);
        right_weight = overlap_transition(right_t, pou_type);

        if i == 1
            weight = right_weight;
        elseif i == Nx
            weight = left_weight;
        else
            weight = min(left_weight, right_weight);
        end
    end

    pu1dx(:, i) = weight;
end

denom = sum(pu1dx, 2);
if any(denom <= eps)
    error('The generated partition of unity has uncovered fine-grid nodes.');
end

pu1dx = pu1dx ./ denom;
end

function y = overlap_transition(t, pou_type)
t = min(1, max(0, t));
switch pou_type
    case 'linear'
        y = t;
    case 'cubic'
        y = t.^2 .* (3 - 2 * t);
    case 'quintic'
        y = t.^3 .* (10 - 15 * t + 6 * t.^2);
    case 'cosine'
        y = 0.5 - 0.5 * cos(pi * t);
    otherwise
        error('Unknown pou_type "%s". Use legacy, linear, cubic, quintic, or cosine.', pou_type);
end
end
