function [applied_field, params] = evaluate_applied_field(t, state_vector, params, step_index)
% EVALUATE_APPLIED_FIELD  Compute applied magnetic field at the current time.
%
%   [applied_field, params] = evaluate_applied_field(t, state_vector, params, step_index)
%
%   Determines the applied boundary magnetic field based on the ramp
%   schedule, computes the physical B-field from the current state (for
%   monitoring), and returns the applied_field struct used by BCs.
%
%   Ramp Schedule (default):
%     - 0 < t <= t_stop * 2/3 : field is ON at full magnitude
%     - t > t_stop * 2/3      : field is OFF
%
%   Inputs:
%     t            : current simulation time
%     state_vector : current state [psi; phi_x; phi_y; phi_z]
%     params       : parameter struct
%     step_index   : integer, current time step (for history storage)
%
%   Outputs:
%     applied_field : struct with fields:
%       .Bx, .By, .Bz : sparse vectors (n_nodes_full × 1)
%                        containing the applied field on boundary nodes.
%     params : updated params (with B-field history stored).
%
%   See also EVALUATE_BFIELD, EVALUATE_RHS

    n_full = params.n_nodes_full;

    % ---------- Compute physical B-field for monitoring ----------
    [Bx_physical, By_physical, Bz_physical] = evaluate_bfield(state_vector, params);
    params.BXT(:, step_index) = Bx_physical;
    params.BYT(:, step_index) = By_physical;
    params.BZT(:, step_index) = Bz_physical;

    % ---------- Determine applied field magnitude ----------
    if t > 0 && t <= params.t_stop * 2/3
        active_Bx = params.applied_Bx;
        active_By = params.applied_By;
        active_Bz = params.applied_Bz;
    else
        active_Bx = 0;
        active_By = 0;
        active_Bz = 0;
    end

    % ---------- Build boundary-node B vectors ----------
    % Boundary nodes where each component is applied
    % Bx is applied on y-faces and z-faces (not x-faces)
    boundary_Bx = [params.y_face_lo_inner, params.z_face_lo_inner, ...
                   params.y_face_hi_inner, params.z_face_hi_inner];
    % By is applied on z-faces and x-faces
    boundary_By = [params.z_face_lo_inner, params.x_face_lo_inner, ...
                   params.z_face_hi_inner, params.x_face_hi_inner];
    % Bz is applied on x-faces and y-faces
    boundary_Bz = [params.x_face_lo_inner, params.y_face_lo_inner, ...
                   params.x_face_hi_inner, params.y_face_hi_inner];

    % Sparse boundary vectors
    Bx_applied = sparse(n_full, 1);
    By_applied = sparse(n_full, 1);
    Bz_applied = sparse(n_full, 1);

    if active_Bx > 0
        params.B_real_x = active_Bx;
        Bx_applied = Bx_applied + sparse(boundary_Bx, 1, active_Bx, n_full, 1);
    else
        params.B_real_x = 0;
    end

    if active_By > 0
        params.B_real_y = active_By;
        By_applied = By_applied + sparse(boundary_By, 1, active_By, n_full, 1);
    else
        params.B_real_y = 0;
    end

    if active_Bz > 0
        params.B_real_z = active_Bz;
        Bz_applied = Bz_applied + sparse(boundary_Bz, 1, active_Bz, n_full, 1);
    else
        params.B_real_z = 0;
    end

    applied_field.Bx = Bx_applied;
    applied_field.By = By_applied;
    applied_field.Bz = Bz_applied;

    params.applied_field = applied_field;
end
