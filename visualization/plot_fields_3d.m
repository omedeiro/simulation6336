function plot_fields_3d(step_index, state_history, params, varargin)
% PLOT_FIELDS_3D  Visualize |ψ|² and B-field components as 3D scatter plots.
%
%   plot_fields_3d(step_index, state_history, params)
%   plot_fields_3d(step_index, state_history, params, 'save_gif', true)
%
%   Produces a four-panel figure (2×2):
%     (1,1) |ψ|² order parameter
%     (1,2) Bx component
%     (2,1) By component
%     (2,2) Bz component
%
%   Optional Name-Value Pairs:
%     'save_gif'   : logical (default false).
%     'gif_dir'    : char (default 'gifs').
%
%   See also EVALUATE_BFIELD, PLOT_FIELDS_2D

    ip = inputParser;
    addParameter(ip, 'save_gif', false, @islogical);
    addParameter(ip, 'gif_dir', 'gifs', @ischar);
    parse(ip, varargin{:});

    Nx = params.Nx;
    Ny = params.Ny;
    Nz = params.Nz;
    hx = params.hx;
    hy = params.hy;
    hz = params.hz;

    n_int  = (Nx - 1) * (Ny - 1) * max(Nz - 1, 1);
    n_int2 = numel(params.bfield_interior);

    % 3D coordinate grids for interior nodes
    [yy, xx, zz] = meshgrid(hx:hx:hx*(Nx-1), hy:hy:hy*(Ny-1), hz:hz:hz*(Nz-1));
    xx = xx(:);  yy = yy(:);  zz = zz(:);

    [yy2, xx2, zz2] = meshgrid(hx:hx:hx*(Nx-2), hy:hy:hy*(Ny-2), hz:hz:hz*(Nz-2));
    xx2 = xx2(:);  yy2 = yy2(:);  zz2 = zz2(:);

    % Compute B-field
    [Bx, By, Bz] = evaluate_bfield(state_history(:, step_index), params);

    % Extract state
    psi   = state_history(1:n_int, step_index);
    psi_sq = abs(psi).^2 ./ max(abs(psi).^2, 1);

    fig = figure(1);
    fig.Position = [100, 200, 1000, 800];
    fig.Color = [1, 1, 1];

    % |ψ|²
    subplot(2, 2, 1)
    scatter3(xx, yy, zz, 36, psi_sq, 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar;  axis equal;  view(0, 90)
    title('|\psi|^2');  xlabel('x');  ylabel('y');  zlabel('z')

    % Bx
    subplot(2, 2, 2)
    bx_norm = real(Bx) ./ max(max(real(Bx)), 1);
    scatter3(xx2, yy2, zz2, 36, bx_norm, 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar;  axis equal;  view(0, 90)
    title('B_x');  xlabel('x');  ylabel('y');  zlabel('z')

    % By
    subplot(2, 2, 3)
    by_norm = real(By) ./ max(max(real(By)), 1);
    scatter3(xx2, yy2, zz2, 36, by_norm, 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar;  axis equal;  view(0, 90)
    title('B_y');  xlabel('x');  ylabel('y');  zlabel('z')

    % Bz
    subplot(2, 2, 4)
    scatter3(xx2, yy2, zz2, 36, real(Bz), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar;  axis equal;  view(0, 90)
    title('B_z');  xlabel('x');  ylabel('y');  zlabel('z')

    sgtitle(sprintf('t = %.4f, applied Bz = %.4f', ...
        params.t(step_index), params.B_real_z_history(min(step_index, end))))
    drawnow

    % Save GIF
    if ip.Results.save_gif
        gif_dir = ip.Results.gif_dir;
        if ~exist(gif_dir, 'dir'), mkdir(gif_dir); end
        if step_index == 1
            filename = fullfile(gif_dir, ['trapGif', datestr(now, 30), '.gif']);
            gif(filename);
            save(strrep(filename, '.gif', '.mat'), 'params');
        else
            gif
        end
    end
end
