function [x, r_norms] = tgcr_matrix_free_trapezoidal(rhs_func, x_state, ...
    params, applied_field, rhs_vector, tol_gcr, max_iter, eps_mf, dt)
% TGCR_MATRIX_FREE_TRAPEZOIDAL  GCR solver for the trapezoidal Jacobian system.
%
%   [x, r_norms] = tgcr_matrix_free_trapezoidal(rhs_func, x_state, ...
%       params, applied_field, rhs_vector, tol_gcr, max_iter, eps_mf, dt)
%
%   Solves the linear system:
%
%       J_trap * x = b
%
%   where J_trap = I - (dt/2) * df/dx is the Jacobian of the trapezoidal
%   residual g(x) = x - (dt/2)*f(x) - gamma.
%
%   The Jacobian-vector product J_trap * p is approximated by:
%
%       J_trap * p ≈ p - (dt/2) * [f(x + ε*p) - f(x)] / ε
%
%   with ε = eps_mf * (1 + ||x||) / ||p||.
%
%   This is the Generalized Conjugate Residual (GCR) method, which is
%   equivalent to GMRES for symmetric systems but works for non-symmetric.
%
%   Inputs:
%     rhs_func      : @(X, p, u) → dX/dt
%     x_state       : current state for Jacobian evaluation
%     params        : parameter struct
%     applied_field : struct with .Bx, .By, .Bz
%     rhs_vector    : right-hand side b of the linear system
%     tol_gcr       : relative residual tolerance
%     max_iter      : maximum GCR iterations
%     eps_mf        : perturbation scale for directional derivative
%     dt            : time step (for trapezoidal J = I - dt/2 * J_f)
%
%   Outputs:
%     x       : solution vector
%     r_norms : vector of ||r_k|| / ||r_0|| at each iteration
%
%   See also NEWTON_GCR_TRAPEZOIDAL, TGCR_MATRIX_FREE

    % Initial guess: x = 0
    x = zeros(size(rhs_vector));
    residual = rhs_vector;
    r_norms(1) = norm(residual, 2);

    k = 0;
    p  = [];  % search directions
    Ap = [];  % J_trap * search directions

    while (r_norms(k + 1) / r_norms(1) > tol_gcr) && (k <= max_iter)
        k = k + 1;

        % New search direction = current residual
        p(:, k) = residual; %#ok<AGROW>

        % ---- Matrix-free J_trap * p(:,k) ----
        epsilon = eps_mf * (1 + norm(x_state)) / norm(p(:, k));
        f_perturbed = feval(rhs_func, x_state + epsilon * p(:, k), params, applied_field);
        f_current   = feval(rhs_func, x_state,                     params, applied_field);
        Ap(:, k) = p(:, k) - (dt / 2) * (f_perturbed - f_current) / epsilon; %#ok<AGROW>

        % ---- Orthogonalize against previous Ap vectors ----
        if k > 1
            for j = 1:(k - 1)
                beta = Ap(:, k)' * Ap(:, j);
                p(:, k)  = p(:, k)  - beta * p(:, j);
                Ap(:, k) = Ap(:, k) - beta * Ap(:, j);
            end
        end

        % ---- Normalize ----
        norm_Ap = norm(Ap(:, k), 2);
        Ap(:, k) = Ap(:, k) / norm_Ap;
        p(:, k)  = p(:, k)  / norm_Ap;

        % ---- Update solution and residual ----
        alpha = residual' * Ap(:, k);
        x = x + alpha * p(:, k);
        residual = residual - alpha * Ap(:, k);

        r_norms(k + 1) = norm(residual, 2); %#ok<AGROW>
    end
end
