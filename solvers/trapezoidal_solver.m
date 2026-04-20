function [state_history, params, final_rhs, errf_k, errDeltax_k] = ...
    trapezoidal_solver(rhs_func, initial_state, params, field_func, ...
                       t_start, t_stop, dt_initial, visualize_func)
% TRAPEZOIDAL_SOLVER  Implicit trapezoidal time integrator with Newton-GCR.
%
%   [state_history, params, final_rhs, errf_k, errDeltax_k] = ...
%       trapezoidal_solver(rhs_func, initial_state, params, field_func, ...
%                          t_start, t_stop, dt_initial, visualize_func)
%
%   Advances the TDGL state using the trapezoidal rule:
%
%       X^{n+1} = X^n + (dt/2) * [f(X^n) + f(X^{n+1})]
%
%   The implicit equation is solved at each step using Newton iteration
%   with a matrix-free GCR (Generalized Conjugate Residual) inner solver.
%
%   If Newton fails to converge, dt is reduced by 10× and the step is
%   retried (adaptive time stepping).
%
%   Inputs:
%     rhs_func       : function handle @(X, p, u) → dX/dt
%     initial_state  : column vector (4*n_interior × 1)
%     params         : struct from setup_parameters + construct_grid_indices
%     field_func     : function handle @(t, X, p, n) → [applied_field, params]
%     t_start        : start time
%     t_stop         : end time
%     dt_initial     : initial time step (may be reduced on non-convergence)
%     visualize_func : function handle @(n, X, p) or [] to skip plotting
%
%   Outputs:
%     state_history : matrix of all accepted states
%     params        : updated parameter struct (with B-field history)
%     final_rhs     : last RHS evaluation
%     errf_k        : last Newton equation-error norm
%     errDeltax_k   : last Newton step-error norm
%
%   See also FORWARD_EULER, NEWTON_GCR_TRAPEZOIDAL, EVALUATE_RHS

    % Newton-GCR tolerances
    tol_gcr_residual  = 1e-4;   % GCR residual convergence
    eps_matrix_free   = 1e-4;   % finite-difference perturbation for Jacobian-vector products
    tol_equation_err  = 1e-3;   % Newton: max |f(x)| for convergence
    tol_step_err      = 1e-3;   % Newton: max |Δx| for convergence
    rel_step_err      = 1;      % relative step tolerance (1 = effectively disabled)
    max_newton_iter   = 20;     % maximum Newton iterations per time step

    % Pre-allocate B-field history
    n_max_steps = ceil((t_stop - t_start) / dt_initial);
    params.BXT = sparse(length(params.bfield_interior), n_max_steps);
    params.BYT = sparse(length(params.bfield_interior), n_max_steps);
    params.BZT = sparse(length(params.bfield_interior), n_max_steps);

    % Initialise
    state_history(:, 1) = initial_state;
    params.t(1) = t_start;
    dt = dt_initial;
    step = 1;

    % ===================== Time Integration Loop =====================
    while params.t(step) < t_stop

        % Evaluate applied field at current time
        [applied_field, params] = feval(field_func, params.t(step), ...
            state_history(:, step), params, step);
        params.B_real_x_history(step) = params.B_real_x;
        params.B_real_y_history(step) = params.B_real_y;
        params.B_real_z_history(step) = params.B_real_z;

        % ---- Explicit predictor (Forward Euler) ----
        f_current = feval(rhs_func, state_history(:, step), params, applied_field);
        x_predicted = state_history(:, step) + dt * f_current;

        % gamma = X^n + (dt/2)*f(X^n)  (the known part of trapezoidal rule)
        gamma = state_history(:, step) + (dt / 2) * f_current;

        % ---- Newton-GCR implicit correction ----
        [x_corrected, converged, errf_k, errDeltax_k] = ...
            newton_gcr_trapezoidal(rhs_func, x_predicted, params, applied_field, ...
                tol_equation_err, tol_step_err, rel_step_err, max_newton_iter, ...
                false, false, [], tol_gcr_residual, eps_matrix_free, gamma, dt);

        if converged
            state_history(:, step + 1) = x_corrected;
            params.t(step + 1) = params.t(step) + dt;
            fprintf('  t = %.4f\n', params.t(step));
            step = step + 1;
        else
            dt = dt / 10;
            fprintf('  Newton did not converge. Reducing dt to %.2e\n', dt);
        end

        % Optional visualization
        if ~isempty(visualize_func)
            feval(visualize_func, step, state_history, params);
        end
    end

    final_rhs = f_current;
end
