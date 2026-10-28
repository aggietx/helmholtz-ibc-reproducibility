function sol = run_fixed_point(A, b, prec, tol, max_pc_calls, max_seconds)
%RUN_FIXED_POINT Run x <- x + B(b-Ax) with true-residual stopping.

if nargin < 6
    max_seconds = inf;
end

x = zeros(size(b));
normb = norm(b);
history = nan(max_pc_calls + 1, 1);
flag = 1;
solve_timer = tic;

for pc_calls = 0:max_pc_calls
    r = b - A * x;
    history(pc_calls + 1) = norm(r) / normb;
    if history(pc_calls + 1) <= tol
        flag = 0;
        break;
    end
    if ~isfinite(history(pc_calls + 1)) || history(pc_calls + 1) > 1e8
        flag = 3;
        break;
    end
    if pc_calls == max_pc_calls
        break;
    end
    if toc(solve_timer) > max_seconds
        flag = 2;
        break;
    end
    x = x + prec(r);
end

sol = struct();
sol.x = x;
sol.flag = flag;
sol.pc_calls = pc_calls;
sol.true_relres = history(pc_calls + 1);
sol.residual_history = history(1:pc_calls + 1);
sol.solve_time = toc(solve_timer);
end
