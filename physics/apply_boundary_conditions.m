function [psi_full, phi_x_full, phi_y_full, phi_z_full] = ...
    apply_boundary_conditions(psi_full_interior, phi_x_full_interior, ...
                              phi_y_full_interior, phi_z_full_interior, ...
                              applied_field, params)
% APPLY_BOUNDARY_CONDITIONS  Apply BCs to expand interior → full-grid fields.
%
%   [psi_full, phi_x_full, phi_y_full, phi_z_full] = ...
%       apply_boundary_conditions(psi_int, phi_x_int, phi_y_int, phi_z_int, ...
%                                 applied_field, params)
%
%   Takes interior-only data (already placed into sparse full-grid vectors
%   at the interior indices) and fills in the boundary nodes using:
%
%     1. Zero-current BCs for psi:
%          psi(face_lo) = psi(first_inner) * exp(-i * phi(face_lo))
%          psi(face_hi) = psi(last_inner)  * exp(+i * phi(last_inner))
%
%     2. Magnetic field BCs for phi components (from applied B × area):
%          phi_y on x-faces:  phi_y(lo) = phi_y(first) - Bz * hx * hy
%                              phi_y(hi) = phi_y(last)  + Bz * hx * hy
%          phi_z on x-faces:  phi_z(lo) = phi_z(first) + By * hz * hx
%                              phi_z(hi) = phi_z(last)  - By * hz * hx
%          (and cyclically for y and z faces)
%
%     3. Periodic BCs (if enabled): boundary = opposite interior copy.
%
%   Inputs:
%     psi_full_interior   : sparse vector (n_nodes_full × 1), interior values set.
%     phi_{x,y,z}_full_interior : same.
%     applied_field       : struct with fields Bx, By, Bz (sparse vectors).
%     params              : struct from setup_parameters + construct_grid_indices.
%
%   Outputs:
%     psi_full, phi_x_full, phi_y_full, phi_z_full :
%         sparse vectors (n_nodes_full × 1) with boundary values filled in.
%
%   See also EVALUATE_RHS, EVALUATE_APPLIED_FIELD

    n_full = params.n_nodes_full;
    hx = params.hx;
    hy = params.hy;
    hz = params.hz;

    % Start from copies with interior values
    psi   = psi_full_interior;
    phi_x = phi_x_full_interior;
    phi_y = phi_y_full_interior;
    phi_z = phi_z_full_interior;

    % Cache "clean" interior copies for BC expressions
    psi0   = psi_full_interior;
    phi_x0 = phi_x_full_interior;
    phi_y0 = phi_y_full_interior;
    phi_z0 = phi_z_full_interior;

    % Set phi at last-cell faces to zero (link variables end at N, not N+1)
    phi_x(params.x_last) = 0;
    phi_y(params.y_last) = 0;
    phi_z(params.z_last) = 0;

    phi_x0(params.x_last) = 0;
    phi_y0(params.y_last) = 0;
    phi_z0(params.z_last) = 0;

    % ===================== X Boundaries =====================
    if params.periodic_x
        % Periodic: face_lo = last_inner, face_hi = first_inner
        psi   = psi   + sparse(params.x_face_lo_inner, 1, psi0(params.x_last_inner),   n_full, 1);
        psi   = psi   + sparse(params.x_face_hi_inner, 1, psi0(params.x_first_inner),  n_full, 1);
        phi_x = phi_x + sparse(params.x_face_lo_inner, 1, phi_x0(params.x_last_inner), n_full, 1);
        phi_x = phi_x + sparse(params.x_face_hi_inner, 1, phi_x0(params.x_first_inner),n_full, 1);
    else
        % Zero-current BCs for psi on x-faces
        psi = psi + sparse(params.x_face_lo_inner, 1, ...
            psi0(params.x_first_inner) .* exp(-1i * phi_x0(params.x_face_lo_inner)), n_full, 1);
        psi = psi + sparse(params.x_face_hi_inner, 1, ...
            psi0(params.x_last_inner) .* exp(1i * phi_x0(params.x_last_inner)), n_full, 1);

        % Magnetic field BCs for phi_y on x-faces
        phi_y = phi_y + sparse(params.x_face_lo_inner, 1, ...
            -applied_field.Bz(params.x_face_lo_inner) * hx * hy + phi_y0(params.x_first_inner), n_full, 1);
        phi_y = phi_y + sparse(params.x_face_hi_inner, 1, ...
            +applied_field.Bz(params.x_face_hi_inner) * hx * hy + phi_y0(params.x_last_inner), n_full, 1);

        % Magnetic field BCs for phi_z on x-faces
        phi_z = phi_z + sparse(params.x_face_lo_inner, 1, ...
            +applied_field.By(params.x_face_lo_inner) * hz * hx + phi_z0(params.x_first_inner), n_full, 1);
        phi_z = phi_z + sparse(params.x_face_hi_inner, 1, ...
            -applied_field.By(params.x_face_hi_inner) * hz * hx + phi_z0(params.x_last_inner), n_full, 1);
    end

    % ===================== Y Boundaries =====================
    if params.periodic_y
        psi   = psi   + sparse(params.y_face_lo_inner, 1, psi0(params.y_last_inner),   n_full, 1);
        psi   = psi   + sparse(params.y_face_hi_inner, 1, psi0(params.y_first_inner),  n_full, 1);
        phi_y = phi_y + sparse(params.y_face_lo_inner, 1, phi_y0(params.y_last_inner), n_full, 1);
        phi_y = phi_y + sparse(params.y_face_hi_inner, 1, phi_y0(params.y_first_inner),n_full, 1);
    else
        % Zero-current BCs for psi on y-faces
        psi = psi + sparse(params.y_face_lo_inner, 1, ...
            psi0(params.y_first_inner) .* exp(-1i * phi_y0(params.y_face_lo_inner)), n_full, 1);
        psi = psi + sparse(params.y_face_hi_inner, 1, ...
            psi0(params.y_last_inner) .* exp(1i * phi_y0(params.y_last_inner)), n_full, 1);

        % Magnetic field BCs for phi_x on y-faces
        phi_x = phi_x + sparse(params.y_face_lo_inner, 1, ...
            +applied_field.Bz(params.y_face_lo_inner) * hx * hy + phi_x0(params.y_first_inner), n_full, 1);
        phi_x = phi_x + sparse(params.y_face_hi_inner, 1, ...
            -applied_field.Bz(params.y_face_hi_inner) * hx * hy + phi_x0(params.y_last_inner), n_full, 1);

        % Magnetic field BCs for phi_z on y-faces
        phi_z = phi_z + sparse(params.y_face_lo_inner, 1, ...
            -applied_field.Bx(params.y_face_lo_inner) * hy * hz + phi_z0(params.y_first_inner), n_full, 1);
        phi_z = phi_z + sparse(params.y_face_hi_inner, 1, ...
            +applied_field.Bx(params.y_face_hi_inner) * hy * hz + phi_z0(params.y_last_inner), n_full, 1);
    end

    % ===================== Z Boundaries =====================
    if params.periodic_z
        psi   = psi   + sparse(params.z_face_lo_inner, 1, psi0(params.z_last_inner),   n_full, 1);
        psi   = psi   + sparse(params.z_face_hi_inner, 1, psi0(params.z_first_inner),  n_full, 1);
        phi_z = phi_z + sparse(params.z_face_lo_inner, 1, phi_y0(params.z_last_inner), n_full, 1);
        phi_z = phi_z + sparse(params.z_face_hi_inner, 1, phi_y0(params.z_first_inner),n_full, 1);
    else
        % Zero-current BCs for psi on z-faces
        psi = psi + sparse(params.z_face_lo_inner, 1, ...
            psi0(params.z_first_inner) .* exp(-1i * phi_z0(params.z_face_lo_inner)), n_full, 1);
        psi = psi + sparse(params.z_face_hi_inner, 1, ...
            psi0(params.z_last_inner) .* exp(1i * phi_z0(params.z_last_inner)), n_full, 1);

        % Magnetic field BCs for phi_x on z-faces
        phi_x = phi_x + sparse(params.z_face_lo_inner, 1, ...
            -applied_field.By(params.z_face_lo_inner) * hz * hx + phi_x0(params.z_first_inner), n_full, 1);
        phi_x = phi_x + sparse(params.z_face_hi_inner, 1, ...
            +applied_field.By(params.z_face_hi_inner) * hz * hx + phi_x0(params.z_last_inner), n_full, 1);

        % Magnetic field BCs for phi_y on z-faces
        phi_y = phi_y + sparse(params.z_face_lo_inner, 1, ...
            +applied_field.Bx(params.z_face_lo_inner) * hy * hz + phi_y0(params.z_first_inner), n_full, 1);
        phi_y = phi_y + sparse(params.z_face_hi_inner, 1, ...
            -applied_field.Bx(params.z_face_hi_inner) * hy * hz + phi_y0(params.z_last_inner), n_full, 1);
    end

    % Return final full-grid fields
    psi_full   = psi;
    phi_x_full = phi_x;
    phi_y_full = phi_y;
    phi_z_full = phi_z;
end
