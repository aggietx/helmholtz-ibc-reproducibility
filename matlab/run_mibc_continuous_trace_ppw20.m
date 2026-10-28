function records = run_mibc_continuous_trace_ppw20(grids)
% Change only the boundary load space; retain the modified IBC local operator.
here = fileparts(mfilename('fullpath'));
addpath(here,'-begin');
if nargin < 1, grids = [64 128 256 512 1024]; end
out = fullfile(here,'data','mibc_continuous_trace_ppw20');
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
    for ell = [0 1]
        F = build_hel_mibc_edge_loads(dims(1),dims(2),ell,Me,1);
        G = build_hel_mibc_edge_loads(dims(1),dims(2),ell,Me,1,'continuous_boundary');
        assert(size(G,2)==4*2^ell && rank(G)==size(G,2));
        assert(norm(F*(F\G)-G,'fro')/norm(G,'fro')<1e-12);
    end
end
fprintf('PASS: load-space inclusion and ranks for square/rectangular patches.\n');
for n = grids
    keys = cellfun(@(r) r.key,records,'UniformOutput',false);
    required = {};
    for ell=[0 1]
        for method={'fgmres20','stationary'}
            required{end+1}=sprintf('%d_%d_%s',n,ell,method{1}); %#ok<AGROW>
        end
    end
    if all(ismember(required,keys)), continue; end
    q=16; k=2*pi*n/20;
    ctx=build_hel_mibc_fine_context(struct('fine_n',n,'coarse_n',n/q, ...
        'k',k,'tao',0,'solve_reference',false,'verbose',1));
    dd=build_hel_dd_info_modified_ibc(ctx,q,struct('verbose',1, ...
        'pou_type','cubic','ibc_coeff',1i/ctx.Hx,'reuse_homogeneous_factors',true));
    for ell=[0 1]
        c=solve_hel_mibc_wavelet_case(ctx,ell,q,struct('verbose',1, ...
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
            records{end+1}=struct('key',key,'n',n,'ell',ell,'k',k,'q',q, ...
                'trace_family','continuous_boundary','ntrace',4*2^ell, ...
                'method',method{1},'solution',s); %#ok<AGROW>
            save(file,'records','-v7.3');
            fprintf('TRACE_RESULT N=%d ell=%d method=%s PC=%d flag=%d residual=%.12e\n', ...
                n,ell,method{1},s.pc_calls,s.flag,s.true_relres);
        end
        clear prec c;
    end
    clear ctx dd;
end
end
