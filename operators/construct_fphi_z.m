function F_phi_z = construct_fphi_z(psi_full, phi_x_full, phi_y_full, phi_z_full, params)
% CONSTRUCT_FPHI_Z  Build the nonlinear forcing term for the phi_z equation.
%
%   F_phi_z = construct_fphi_z(psi_full, phi_x_full, phi_y_full, phi_z_full, params)
%
%   Computes the curl-curl and supercurrent contributions to dphi_z/dt:
%
%       F_phi_z(m) = kappa^2/hx^2 * [curly contribution from phi_x]
%                  + kappa^2/hy^2 * [curlx contribution from phi_y]
%                  + Im[exp(-i*phi_z(m)) * conj(psi(m)) * psi(m+stride_k)]
%
%   Returns zero for 2D simulations (Nz == 1).
%
%   Inputs:
%     psi_full, phi_x_full, phi_y_full, phi_z_full :
%         Field arrays on the full grid (length n_nodes_full).
%     params : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     F_phi_z : sparse vector (n_interior × 1).
%
%   See also CONSTRUCT_FPHI_X, CONSTRUCT_FPHI_Y, CONSTRUCT_FPSI

    m_compact = params.interior_numbering;

    if params.Nz <= 1
        % 2D case: no z-dynamics
        F_phi_z = sparse(m_compact, 1, 0);
        return
    end

    stride_j  = params.stride_j;
    stride_k  = params.stride_k;
    kappa     = params.kappa;
    hx        = params.hx;
    hy        = params.hy;
    m_full    = params.interior_to_full;

    % Supercurrent: Im[exp(-i*phi_z(m)) * psi*(m) * psi(m+stride_k)]
    supercurrent = imag(exp(-1i * phi_z_full(m_full)) .* ...
                        conj(psi_full(m_full)) .* psi_full(m_full + stride_k));

    % Curl-curl from phi_x
    curl_x = (kappa^2 / hx^2) * ( ...
        -phi_x_full(m_full + stride_k) + phi_x_full(m_full) + ...
         phi_x_full(m_full + stride_k - 1) - phi_x_full(m_full - 1));

    % Curl-curl from phi_y
    curl_y = (kappa^2 / hy^2) * ( ...
        -phi_y_full(m_full + stride_k) + phi_y_full(m_full) + ...
         phi_y_full(m_full + stride_k - stride_j) - phi_y_full(m_full - stride_j));

    values = curl_x + curl_y + supercurrent;

    F_phi_z = sparse(m_compact, 1, values);
end
