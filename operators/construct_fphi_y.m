function F_phi_y = construct_fphi_y(psi_full, phi_x_full, phi_y_full, phi_z_full, params)
% CONSTRUCT_FPHI_Y  Build the nonlinear forcing term for the phi_y equation.
%
%   F_phi_y = construct_fphi_y(psi_full, phi_x_full, phi_y_full, phi_z_full, params)
%
%   Computes the curl-curl and supercurrent contributions to dphi_y/dt:
%
%       F_phi_y(m) = kappa^2/hz^2 * [curlx contribution from phi_z]  (3D only)
%                  + kappa^2/hx^2 * [curlz contribution from phi_x]
%                  + Im[exp(-i*phi_y(m)) * conj(psi(m)) * psi(m+stride_j)]
%
%   Inputs:
%     psi_full, phi_x_full, phi_y_full, phi_z_full :
%         Field arrays on the full grid (length n_nodes_full).
%     params : struct from setup_parameters + construct_grid_indices.
%
%   Output:
%     F_phi_y : sparse vector (n_interior × 1).
%
%   See also CONSTRUCT_FPHI_X, CONSTRUCT_FPHI_Z, CONSTRUCT_FPSI

    stride_j  = params.stride_j;
    stride_k  = params.stride_k;
    kappa     = params.kappa;
    hx        = params.hx;
    hz        = params.hz;
    m_full    = params.interior_to_full;
    m_compact = params.interior_numbering;

    % Supercurrent: Im[exp(-i*phi_y(m)) * psi*(m) * psi(m+stride_j)]
    supercurrent = imag(exp(-1i * phi_y_full(m_full)) .* ...
                        conj(psi_full(m_full)) .* psi_full(m_full + stride_j));

    % Curl-curl from phi_x (always present)
    curl_x = (kappa^2 / hx^2) * ( ...
        -phi_x_full(m_full + stride_j) + phi_x_full(m_full) + ...
         phi_x_full(m_full + stride_j - 1) - phi_x_full(m_full - 1));

    if params.Nz > 1
        % Curl-curl from phi_z
        curl_z = (kappa^2 / hz^2) * ( ...
            -phi_z_full(m_full + stride_j) + phi_z_full(m_full) + ...
             phi_z_full(m_full + stride_j - stride_k) - phi_z_full(m_full - stride_k));

        values = curl_x + curl_z + supercurrent;
    else
        values = curl_x + supercurrent;
    end

    F_phi_y = sparse(m_compact, 1, values);
end
