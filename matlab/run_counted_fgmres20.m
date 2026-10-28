function sol = run_counted_fgmres20(A, b, tol, max_pc_calls, max_seconds, prec)
%RUN_COUNTED_FGMRES20 Right-preconditioned FGMRES(20) with true-residual stopping.

restart = 20;
bnorm = norm(b);
if bnorm == 0
    sol = struct('x', zeros(size(b)), 'flag', 0, 'relres', 0, ...
        'iter', 0, 'resvec', 0, 'pc_calls', 0, 'solve_time', 0, ...
        'true_rel_res', 0);
    return;
end

x = zeros(size(b));
pc_calls = 0;
resvec = norm(b - A * x) / bnorm;
solve_timer = tic;
flag = 1;

while pc_calls < max_pc_calls
    r = b - A * x;
    beta = norm(r);
    true_rel_res = beta / bnorm;
    if true_rel_res <= tol
        flag = 0;
        break;
    end

    cycle_length = min(restart, max_pc_calls - pc_calls);
    V = zeros(numel(b), cycle_length + 1, 'like', b);
    Z = zeros(numel(b), cycle_length, 'like', b);
    H = zeros(cycle_length + 1, cycle_length, 'like', b);
    g = zeros(cycle_length + 1, 1, 'like', b);
    V(:, 1) = r / beta;
    g(1) = beta;

    for j = 1:cycle_length
        if toc(solve_timer) > max_seconds
            error('hel_mibc:solveTimeLimit', ...
                'FGMRES(20) exceeded the %.1f second solve limit.', max_seconds);
        end

        Z(:, j) = prec(V(:, j));
        pc_calls = pc_calls + 1;
        w = A * Z(:, j);

        % Classical Gram--Schmidt matches the PETSc FGMRES configuration.
        H(1:j, j) = V(:, 1:j)' * w;
        w = w - V(:, 1:j) * H(1:j, j);
        hnext = norm(w);
        H(j + 1, j) = hnext;
        if hnext > eps(class(H)) * norm(H(1:j, j))
            V(:, j + 1) = w / hnext;
        end

        y = H(1:j + 1, 1:j) \ g(1:j + 1);
        x_trial = x + Z(:, 1:j) * y;
        true_rel_res = norm(A * x_trial - b) / bnorm;
        resvec(end + 1, 1) = true_rel_res; %#ok<AGROW>
        if true_rel_res <= tol
            x = x_trial;
            flag = 0;
            break;
        end

        if hnext == 0
            x = x_trial;
            break;
        end
    end

    if flag == 0
        break;
    end
    x = x_trial;
end

true_rel_res = norm(A * x - b) / bnorm;
sol = struct('x', x, 'flag', flag, 'relres', true_rel_res, ...
    'iter', pc_calls, 'resvec', resvec, 'pc_calls', pc_calls, ...
    'solve_time', toc(solve_timer), 'true_rel_res', true_rel_res);
end
