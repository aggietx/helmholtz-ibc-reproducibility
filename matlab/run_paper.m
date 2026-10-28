function run_paper(experiment, grids)
%RUN_PAPER Reproduce one manuscript table; checkpoints are saved in matlab/data.
if nargin < 2, grids = [64 128 256 512]; end
addpath(fileparts(mfilename('fullpath')));
switch experiment
    case 'asymptotic'
        run_mibc_half_hthreehalves(grids,0.01,2.01,0.5);
    case 'hthreehalves'
        run_mibc_half_hthreehalves(grids,0.01,1.5,0.5);
    case 'ppw20'
        run_mibc_continuous_trace_ppw20(grids);
    otherwise
        error('Choose asymptotic, hthreehalves, or ppw20.');
end
end
