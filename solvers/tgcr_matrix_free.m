function [x, r_norms] = tgcr_matrix_free(rhs_func, x_state, params, ...
    applied_field, rhs_vector, tol_gcr, max_iter, eps_mf)
% TGCR_MATRIX_FREE  GCR solver for the Jacobian system df/dx * x = b.
%
%   [x, r_norms] = tgcr_matrix_free(rhs_func, x_state, params, ...
%       applied_field, rhs_vector, tol_gcr, max_iter, eps_mf)
%
%   Solves the linear system:
%
%       (df/dx) * x = b
%
%   where df/dx is the Jacobian of the TDGL RHS f. The Jacobian-vector
%   product is approximated via finite differences:
%
%       (df/dx) * p ≈ [f(x + ε*p) - f(x)] / ε
%
%   with ε = eps_mf * (1 + ||x||) / ||p||.
%
%   This solver is used by Newton methods for steady-state problems or
%   non-trapezoidal implicit methods.
%
%   Inputs:
%     rhs_func      : @(X, p, u) → dX/dt
%     x_state       : current state for Jacobian evaluation
%     params        : parameter struct
%     applied_field : struct with .Bx, .By, .Bz
%     rhs_vector    : right-hand side b
%     tol_gcr       : relative residual tolerance
%     max_iter      : maximum GCR iterations
%     eps_mf        : perturbation scale for directional derivative
%
%   Outputs:
%     x       : solution vector
%     r_norms : normalized residual norms at each iteration
%
%   See also TGCR_MATRIX_FREE_TRAPEZOIDAL

    x = zeros(size(rhs_vector));
    residual = rhs_vector;
    r_norms(1) = norm(residual, 2);

    k = 0;
    p  = [];
    Ap = [];

    while (r_norms(k + 1) / r_norms(1) > tol_gcr) && (k <= max_iter)
        k = k + 1;
        p(:, k) = residual; %#ok<AGROW>

        % Matrix-free Jacobian-vector product
        epsilon = eps_mf * (1 + norm(x_state)) / norm(p(:, k));
        f_perturbed = feval(rhs_func, x_state + epsilon * p(:, k), params, applied_field);
        f_current   = feval(rhs_func, x_state,                     params, applied_field);
        Ap(:, k) = (f_perturbed - f_current) / epsilon; %#ok<AGROW>

        % Orthogonalize
        if k > 1
            for j = 1:(k - 1)
                beta = Ap(:, k)' * Ap(:, j);
                p(:, k)  = p(:, k)  - beta * p(:, j);
                Ap(:, k) = Ap(:, k) - beta * Ap(:, j);
            end
        end

        % Normalize
        norm_Ap = norm(Ap(:, k), 2);
        Ap(:, k) = Ap(:, k) / norm_Ap;
        p(:, k)  = p(:, k)  / norm_Ap;

        % Update
        alpha = residual' * Ap(:, k);
        x = x + alpha * p(:, k);
        residual = residual - alpha * Ap(:, k);

        r_norms(k + 1) = norm(residual, 2); %#ok<AGROW>
    end

    if r_norms(k + 1) > (tol_gcr * r_norms(1))
        fprintf('  GCR did NOT converge after %d iterations.\n', k);
        x = [];
    end

    r_norms = r_norms / r_norms(1);
end
