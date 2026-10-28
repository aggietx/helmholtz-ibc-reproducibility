function records = run_mibc_half_hthreehalves(grids,epsilon,fine_exponent,overlap_ratio)
% Change only the boundary load space; retain the modified IBC local operator.
here = fileparts(mfilename('fullpath'));
addpath(here,'-begin');
if nargin < 1, grids = [64 128 256 512]; end
if nargin < 2, epsilon = 0.01; end
assert(isscalar(epsilon) && epsilon>0 && epsilon<=0.25);
if nargin < 3, fine_exponent=1.5; end
assert(isscalar(fine_exponent) && fine_exponent>0);
if nargin < 4, overlap_ratio=0.5; end
assert(isscalar(overlap_ratio) && overlap_ratio>0 && overlap_ratio<=1);
out = fullfile(here,'data','mibc_continuous_trace_asymptotic');
if epsilon~=0.25
    out=[out '_epsilon_' strrep(sprintf('%.8g',epsilon),'.','p')];
end
if fine_exponent~=2+epsilon
    out=[out '_h_power_' strrep(sprintf('%.8g',fine_exponent),'.','p')];
end
if overlap_ratio~=1
    out=[out '_deltaH_' strrep(sprintf('%.8g',overlap_ratio),'.','p')];
end
out=[out '_atomic_v7'];
if ~exist(out,'dir'), mkdir(out); end
diary(fullfile(out,'run.log'));
cleanup = onCleanup(@() diary('off')); %#ok<NASGU>
file = fullfile(out,'results.mat');
records = {};
if exist(file,'file'), saved=load(file,'records'); records=saved.records; end
% Check the new load space is a subspace of the independent-edge space,
% including rectangular patches encountered at the physical boundary.
for dims = [32 48; 48 32; 48 48]'
    [~,~,~,~,Me] = get_element_blend_matrix(1/64,1/64,1,0);
    for ell = [0 1 2 3]
        F = build_hel_mibc_edge_loads(dims(1),dims(2),ell,Me,1);
        G = build_hel_mibc_edge_loads(dims(1),dims(2),ell,Me,1,'continuous_boundary');
        assert(size(G,2)==4*2^ell && rank(G)==size(G,2));
        assert(norm(F*(F\G)-G,'fro')/norm(G,'fro')<1e-12);
    end
end
fprintf('PASS: load-space inclusion and ranks for square/rectangular patches.\n');
for n = grids
    k=n^(1/fine_exponent); levels=ceil(log2(k^(0.5+epsilon)));
    q=max(2^ceil(log2(n*k^(-1-epsilon))),2^levels);
    over=overlap_ratio*q;
    assert(over==round(over) && mod(n,q)==0, 'Nonintegral overlap or coarse grid.');
    if overlap_ratio==0.5 && fine_exponent==2+epsilon
        assert(q==2^ceil(log2(n*k^(-1-epsilon))), 'Screenshot coarse-grid rule differs.');
    end
    max_level=log2(q+2*over);
    fprintf('PARAMETERS N=%d epsilon=%.8g k=%.12g ppw=%.8g H/h=%d ell=%d\n', ...
        n,epsilon,k,2*pi*n/k,q,levels);
    fprintf('OVERLAP delta/H=%.8g delta/h=%d interior_width/h=%d L=%.8g\n', ...
        overlap_ratio,over,q+2*over,max_level);
    keys = cellfun(@(r) r.key,records,'UniformOutput',false);
    required = {};
    for ell=levels
        for method={'fgmres20','stationary'}
            required{end+1}=sprintf('%d_%d_%s',n,ell,method{1}); %#ok<AGROW>
        end
    end
    if all(ismember(required,keys)), continue; end
    ctx=build_hel_mibc_fine_context(struct('fine_n',n,'coarse_n',n/q, ...
        'k',k,'tao',0,'solve_reference',false,'verbose',1));
    dd=build_hel_dd_info_modified_ibc(ctx,over,struct('verbose',1, ...
        'pou_type','cubic','ibc_coeff',1i/ctx.Hx,'reuse_homogeneous_factors',true));
    for ell=levels
        c=solve_hel_mibc_wavelet_case(ctx,ell,over,struct('verbose',1, ...
            'trace_family','continuous_boundary','pou_type','cubic', ...
            'keep_basis',1,'keep_coarse_matrix',1,'keep_coarse_factor',1,'dd_info',dd));
        prec=@(r) dd_pre_hel_mibc_2level_hybrid_batched( ...
            c.phiK,c.Phi,c.Uc,c.Lc,c.Pc,c.Qc,dd,r);
        for method={'fgmres20','stationary'}
            key=sprintf('%d_%d_%s',n,ell,method{1});
            if any(strcmp(keys,key)), continue; end
            if strcmp(method{1},'fgmres20')
                s=run_counted_fgmres20(ctx.Astiffd,ctx.Fd,1e-6,200,600,prec);
                s.true_relres=s.true_rel_res;
            else
                s=run_fixed_point(ctx.Astiffd,ctx.Fd,prec,1e-6,200,600);
            end
            s=rmfield(s,'x');
            records{end+1}=struct('key',key,'n',n,'ell',ell,'k',k,'q',q,'epsilon',epsilon, ...
                'fine_exponent',fine_exponent, ...
                'overlap_ratio',overlap_ratio,'overlap_elements',over, ...
                'interior_max_level',max_level, ...
                'trace_family','continuous_boundary','ntrace',4*2^ell, ...
                'method',method{1},'solution',s); %#ok<AGROW>
            checkpoint=[tempname(out) '.mat'];
            save(checkpoint,'records','-v7');
            check=load(checkpoint,'records');
            assert(numel(check.records)==numel(records));
            [ok,msg]=movefile(checkpoint,file,'f');
            assert(ok,msg);
            fprintf('TRACE_RESULT N=%d ell=%d method=%s PC=%d flag=%d residual=%.12e\n', ...
                n,ell,method{1},s.pc_calls,s.flag,s.true_relres);
        end
        clear prec c;
    end
    clear ctx dd;
end
end
