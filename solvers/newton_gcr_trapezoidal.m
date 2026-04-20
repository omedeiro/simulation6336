function [x, converged, errf_k, errDeltax_k, relDeltax_k, iterations] = ...
    newton_gcr_trapezoidal(rhs_func, x0, params, applied_field, ...
        tol_equation, tol_step, tol_relative, max_iter, ...
        visualize, use_finite_diff, jacobian_func, ...
        tol_gcr, eps_mf, gamma, dt)
% NEWTON_GCR_TRAPEZOIDAL  Newton solver with matrix-free GCR for trapezoidal rule.
%
%   [x, converged, errf_k, errDeltax_k, relDeltax_k, iterations] = ...
%       newton_gcr_trapezoidal(rhs_func, x0, params, applied_field, ...
%           tol_equation, tol_step, tol_relative, max_iter, ...
%           visualize, use_finite_diff, jacobian_func, ...
%           tol_gcr, eps_mf, gamma, dt)
%
%   Solves the implicit trapezoidal equation:
%
%       g(x) = x - (dt/2)*f(x) - gamma = 0
%
%   where gamma = X^n + (dt/2)*f(X^n), using Newton iteration with the
%   GCR linear solver (matrix-free Jacobian-vector products via finite
%   differences).
%
%   Convergence is declared when ALL of:
%     - |g(x)| < tol_equation
%     - |Δx|   < tol_step
%     - |Δx|/max|x| < tol_relative
%
%   Inputs:
%     rhs_func        : @(X, p, u) → dX/dt
%     x0              : initial guess (usually explicit predictor)
%     params          : parameter struct
%     applied_field   : struct with .Bx, .By, .Bz
%     tol_equation    : absolute equation error tolerance
%     tol_step        : absolute step error tolerance
%     tol_relative    : relative step error tolerance
%     max_iter        : maximum Newton iterations
%     visualize       : logical, show intermediate results
%     use_finite_diff : logical, use FD Jacobian instead of given
%     jacobian_func   : handle for explicit Jacobian (or [])
%     tol_gcr         : GCR residual tolerance
%     eps_mf          : perturbation for matrix-free directional derivative
%     gamma           : trapezoidal known vector = X^n + (dt/2)*f(X^n)
%     dt              : time step
%
%   Outputs:
%     x, converged, errf_k, errDeltax_k, relDeltax_k, iterations
%
%   See also TRAPEZOIDAL_SOLVER, TGCR_MATRIX_FREE_TRAPEZOIDAL

    N = length(x0);
    max_gcr_iter = max(N, round(N * 0.2));

    % Initialize
    k = 1;
    X_history = x0;
    f = feval(rhs_func, X_history(:, k), params, applied_field);
    errf_k = norm(f, inf);
    errDeltax_k = 0;
    relDeltax_k = 0;

    while k <= max_iter && ...
          (errf_k > tol_equation || errDeltax_k > tol_step || relDeltax_k > tol_relative)

        f = feval(rhs_func, X_history(:, k), params, applied_field);

        % Trapezoidal residual: g(x) = x - (dt/2)*f(x) - gamma
        g_trap = X_history(:, k) - (dt / 2) * f - gamma;

        if ~isempty(eps_mf)
            % Matrix-free GCR solve for Δx in: J_trap * Δx = -g_trap
            delta_x = tgcr_matrix_free_trapezoidal( ...
                rhs_func, X_history(:, k), params, applied_field, ...
                -g_trap, tol_gcr, max_gcr_iter, eps_mf, dt);
        else
            if use_finite_diff
                error('newton_gcr_trapezoidal:notImplemented', ...
                    'Explicit FD Jacobian not supported; use eps_mf instead.');
            else
                Jf = feval(jacobian_func, X_history(:, k), params, applied_field);
                delta_x = tgcr(Jf, -f, tol_gcr, max_gcr_iter);
            end
        end

        % Update
        X_history(:, k + 1) = X_history(:, k) + delta_x;
        k = k + 1;

        % Re-evaluate for convergence check
        f = feval(rhs_func, X_history(:, k), params, applied_field);
        g_trap = X_history(:, k) - (dt / 2) * f - gamma;

        errf_k      = norm(g_trap, inf);
        errDeltax_k = norm(delta_x, inf);
        relDeltax_k = norm(delta_x, inf) / max(abs(X_history(:, k)));
    end

    x = X_history(:, k);
    iterations = k - 1;

    converged = (errf_k <= tol_equation) && ...
                (errDeltax_k <= tol_step) && ...
                (relDeltax_k <= tol_relative);

    if ~converged
        fprintf('  Newton did NOT converge after %d iterations.\n', iterations);
    end
end
