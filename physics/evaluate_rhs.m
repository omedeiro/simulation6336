function F = evaluate_rhs(state_vector, params, applied_field)
% EVALUATE_RHS  Compute the right-hand side dX/dt of the TDGL equations.
%
%   F = evaluate_rhs(state_vector, params, applied_field)
%
%   Evaluates the full RHS of the discretised TDGL system:
%
%       dpsi/dt   = (1/hx^2 L_psi_x + 1/hy^2 L_psi_y + 1/hz^2 L_psi_z) * psi
%                   + (1 - |psi|^2) * psi
%
%       dphi_x/dt = (L_phi_y + L_phi_z) * phi_x + F_phi_x
%       dphi_y/dt = (L_phi_x + L_phi_z) * phi_y + F_phi_y
%       dphi_z/dt = (L_phi_x + L_phi_y) * phi_z + F_phi_z
%
%   Steps:
%     1. Unpack state vector into [psi; phi_x; phi_y; phi_z] interior arrays.
%     2. Place interior values onto full grid via sparse indexing.
%     3. Apply boundary conditions (zero-current + magnetic field).
%     4. Build all sparse operators (LPSI, LPHI, FPSI, FPHI).
%     5. Strip boundary rows from operators.
%     6. Assemble F = [dpsi/dt; dphi_x/dt; dphi_y/dt; dphi_z/dt].
%
%   Inputs:
%     state_vector  : column vector of length 4*n_interior =
%                     [psi_int; phi_x_int; phi_y_int; phi_z_int].
%     params        : struct from setup_parameters + construct_grid_indices.
%     applied_field : struct with .Bx, .By, .Bz (sparse full-grid vectors).
%
%   Output:
%     F : column vector (4*n_interior × 1), time derivative of state.
%
%   See also APPLY_BOUNDARY_CONDITIONS, SETUP_PARAMETERS

    n_int  = params.n_interior;
    n_full = params.n_nodes_full;
    m      = params.interior_to_full;

    % ==================== 1. Unpack State ====================
    psi_int   = state_vector(1         : n_int);
    phi_x_int = state_vector(n_int + 1 : 2 * n_int);
    phi_y_int = state_vector(2*n_int+1 : 3 * n_int);
    phi_z_int = state_vector(3*n_int+1 : 4 * n_int);

    % ==================== 2. Place on Full Grid ====================
    psi_full   = sparse(m, 1, psi_int,   n_full, 1);
    phi_x_full = sparse(m, 1, phi_x_int, n_full, 1);
    phi_y_full = sparse(m, 1, phi_y_int, n_full, 1);
    phi_z_full = sparse(m, 1, phi_z_int, n_full, 1);

    % ==================== 3. Apply Boundary Conditions ====================
    [psi_full, phi_x_full, phi_y_full, phi_z_full] = ...
        apply_boundary_conditions(psi_full, phi_x_full, phi_y_full, phi_z_full, ...
                                  applied_field, params);

    % ==================== 4. Build Operators ====================
    % Covariant Laplacian for psi (state-dependent via Peierls phases)
    L_psi_x = construct_lpsi_x(phi_x_full, params);
    L_psi_y = construct_lpsi_y(phi_y_full, params);
    L_psi_z = construct_lpsi_z(phi_z_full, params);

    % Laplacians for phi (constant, use precomputed if available)
    if isfield(params, 'L_phi_x_precomputed')
        L_phi_x = params.L_phi_x_precomputed;
        L_phi_y = params.L_phi_y_precomputed;
        L_phi_z = params.L_phi_z_precomputed;
    else
        L_phi_x = construct_lphi_x(params);
        L_phi_y = construct_lphi_y(params);
        L_phi_z = construct_lphi_z(params);
    end

    % Nonlinear forcing terms
    F_psi   = construct_fpsi(psi_full, params);
    F_phi_x = construct_fphi_x(psi_full, phi_x_full, phi_y_full, phi_z_full, params);
    F_phi_y = construct_fphi_y(psi_full, phi_x_full, phi_y_full, phi_z_full, params);
    F_phi_z = construct_fphi_z(psi_full, phi_x_full, phi_y_full, phi_z_full, params);

    % ==================== 5. Strip Boundary Rows ====================
    % Remove all-zero rows (boundary nodes have no equations)
    L_psi_x(~any(L_psi_x, 2), :) = [];
    L_psi_y(~any(L_psi_y, 2), :) = [];
    L_psi_z(~any(L_psi_z, 2), :) = [];

    L_phi_x(~any(L_phi_x, 2), :) = [];
    L_phi_y(~any(L_phi_y, 2), :) = [];
    L_phi_z(~any(L_phi_z, 2), :) = [];

    % ==================== 6. Assemble RHS ====================
    hx = params.hx;
    hy = params.hy;
    hz = params.hz;

    dpsi_dt   = (L_psi_x / hx^2 + L_psi_y / hy^2 + L_psi_z / hz^2) * psi_full + F_psi;
    dphi_x_dt = (L_phi_y + L_phi_z) * phi_x_full + F_phi_x;
    dphi_y_dt = (L_phi_x + L_phi_z) * phi_y_full + F_phi_y;
    dphi_z_dt = (L_phi_x + L_phi_y) * phi_z_full + F_phi_z;

    F = [dpsi_dt; dphi_x_dt; dphi_y_dt; dphi_z_dt];
end
