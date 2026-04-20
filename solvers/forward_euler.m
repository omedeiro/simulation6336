function state_history = forward_euler(rhs_func, initial_state, params, ...
                                       field_func, t_start, t_stop, dt, ...
                                       visualize_func)
% FORWARD_EULER  Explicit Forward Euler time integrator for TDGL.
%
%   state_history = forward_euler(rhs_func, initial_state, params, ...
%                                  field_func, t_start, t_stop, dt, ...
%                                  visualize_func)
%
%   Advances the TDGL state from t_start to t_stop using explicit Euler:
%
%       X(:, n+1) = X(:, n) + dt * f(X(:, n))
%
%   where f is the RHS of the TDGL equations.
%
%   Inputs:
%     rhs_func       : function handle @(X, p, u) → dX/dt
%     initial_state  : column vector (4*n_interior × 1)
%     params         : struct from setup_parameters + construct_grid_indices
%     field_func     : function handle @(t, X, p, n) → [applied_field, params]
%     t_start        : start time
%     t_stop         : end time
%     dt             : time step
%     visualize_func : function handle @(n, X, p) or [] to skip plotting
%
%   Output:
%     state_history : matrix (4*n_interior × n_steps+1), all states.
%
%   See also TRAPEZOIDAL_SOLVER, EVALUATE_RHS, EVALUATE_APPLIED_FIELD

    n_steps = ceil((t_stop - t_start) / dt);

    % Pre-allocate history
    state_history = zeros(length(initial_state), n_steps + 1);
    state_history(:, 1) = initial_state;

    time = zeros(1, n_steps + 1);
    time(1) = t_start;

    for n = 1:n_steps
        % Adaptive final step
        dt_actual = min(dt, t_stop - time(n));
        time(n + 1) = time(n) + dt_actual;

        % Evaluate applied field
        [applied_field, params] = feval(field_func, time(n), ...
            state_history(:, n), params, n);

        % Store applied field magnitude history
        params.B_real_x_history(n) = params.B_real_x;
        params.B_real_y_history(n) = params.B_real_y;
        params.B_real_z_history(n) = params.B_real_z;

        % Compute RHS
        f = feval(rhs_func, state_history(:, n), params, applied_field);

        % Forward Euler step
        state_history(:, n + 1) = state_history(:, n) + dt_actual * f;

        % Optional visualization
        if ~isempty(visualize_func)
            feval(visualize_func, n, state_history, params);
        end
    end
end
