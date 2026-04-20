function [Bx, By, Bz] = evaluate_bfield(state_matrix, params)
% EVALUATE_BFIELD  Compute the magnetic field B = curl(A) from the state.
%
%   [Bx, By, Bz] = evaluate_bfield(state_matrix, params)
%
%   Computes the discrete curl of the vector potential [phi_x, phi_y, phi_z]
%   at interior nodes (one layer further inward than the state interior):
%
%       Bx = 1/(hy*hz) * (phi_y(m) - phi_y(m+stride_j_int) - phi_z(m) + phi_z(m+stride_i_int))
%       By = 1/(hz*hx) * (phi_z(m) - phi_z(m+1)            - phi_x(m) + phi_x(m+stride_j_int))
%       Bz = 1/(hx*hy) * (phi_x(m) - phi_x(m+stride_i_int) - phi_y(m) + phi_y(m+1))
%
%   where the strides are computed in the compact interior numbering.
%
%   Inputs:
%     state_matrix : matrix (4*n_interior × n_times). Each column is a
%                    state vector. If n_times > 1, B is computed for each.
%     params       : struct from setup_parameters + construct_grid_indices.
%
%   Outputs:
%     Bx, By, Bz : matrices (n_bfield_interior × n_times).
%
%   See also EVALUATE_APPLIED_FIELD, EVALUATE_RHS

    Nx = params.Nx;
    Ny = params.Ny;
    Nz = params.Nz;
    hx = params.hx;
    hy = params.hy;
    hz = params.hz;

    n_times = size(state_matrix, 2);

    if Nz > 1
        n_int    = (Nx - 1) * (Ny - 1) * (Nz - 1);
        n_bfield = numel(params.bfield_interior);
        stride_j_int = Nx - 1;     % j-stride in compact interior numbering
        stride_k_int = (Nx - 1) * (Ny - 1);  % k-stride

        Bx = sparse(n_bfield, n_times);
        By = sparse(n_bfield, n_times);
        Bz = sparse(n_bfield, n_times);
        m_b = params.bfield_interior;

        for t_idx = 1:n_times
            state = state_matrix(:, t_idx);
            phi_x = state(n_int + 1 : 2 * n_int);
            phi_y = state(2*n_int+1 : 3 * n_int);
            phi_z = state(3*n_int+1 : 4 * n_int);

            Bx(:, t_idx) = (1 / (hy * hz)) * ...
                (phi_y(m_b) - phi_y(m_b + stride_k_int) - ...
                 phi_z(m_b) + phi_z(m_b + stride_j_int));

            By(:, t_idx) = (1 / (hz * hx)) * ...
                (phi_z(m_b) - phi_z(m_b + 1) - ...
                 phi_x(m_b) + phi_x(m_b + stride_k_int));

            Bz(:, t_idx) = (1 / (hx * hy)) * ...
                (phi_x(m_b) - phi_x(m_b + stride_j_int) - ...
                 phi_y(m_b) + phi_y(m_b + 1));
        end
    else
        n_int    = (Nx - 1) * (Ny - 1);
        n_bfield = numel(params.bfield_interior);
        stride_j_int = Nx - 1;

        Bx = sparse(n_bfield, n_times);
        By = sparse(n_bfield, n_times);
        Bz = sparse(n_bfield, n_times);
        m_b = params.bfield_interior;

        for t_idx = 1:n_times
            state = state_matrix(:, t_idx);
            phi_x = state(n_int + 1 : 2 * n_int);
            phi_y = state(2*n_int+1 : 3 * n_int);

            Bx(:, t_idx) = 0;
            By(:, t_idx) = 0;
            Bz(:, t_idx) = (1 / (hx * hy)) * ...
                (phi_x(m_b) - phi_x(m_b + stride_j_int) - ...
                 phi_y(m_b) + phi_y(m_b + 1));
        end
    end
end
