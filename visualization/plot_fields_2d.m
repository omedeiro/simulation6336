function plot_fields_2d(state_history, params, varargin)
% PLOT_FIELDS_2D  Visualize |ψ|² and Bz as 2D image plots across time.
%
%   plot_fields_2d(state_history, params)
%   plot_fields_2d(state_history, params, 'save_gif', true)
%
%   Produces a two-panel figure at evenly spaced time frames:
%     Top:    |ψ|² order parameter (normalized to [0,1])
%     Bottom: Bz magnetic field component
%
%   For 3D simulations, displays a single z-slice (set via params.slice_z).
%
%   Optional Name-Value Pairs:
%     'save_gif'   : logical (default false). Save frames as animated GIF.
%     'gif_dir'    : char (default 'gifs'). Directory for GIF output.
%     'bz_clim'    : [lo hi] (default [-0.03 0.03]). Color limits for Bz.
%     'n_frames'   : integer (default params.frames). Number of frames.
%
%   See also EVALUATE_BFIELD, PLOT_FIELDS_3D

    ip = inputParser;
    addParameter(ip, 'save_gif', false, @islogical);
    addParameter(ip, 'gif_dir', 'gifs', @ischar);
    addParameter(ip, 'bz_clim', [-0.03 0.03], @(x) numel(x)==2);
    addParameter(ip, 'n_frames', params.frames, @isnumeric);
    parse(ip, varargin{:});

    Nx = params.Nx;
    Ny = params.Ny;
    Nz = params.Nz;
    hx = params.hx;
    hy = params.hy;

    n_int  = (Nx - 1) * (Ny - 1) * max(Nz - 1, 1);
    n_int2 = (Nx - 2) * (Ny - 2) * max(Nz - 2, 1);
    mk     = (Nx - 1) * (Ny - 1);
    mk2    = (Nx - 2) * (Ny - 2);

    % Coordinate grids for plotting
    x_coords = hx : hx : hx * (Nx - 1);
    y_coords = hy : hy : hy * (Ny - 1);
    x_coords2 = hx : hx : hx * (Nx - 2);
    y_coords2 = hy : hy : hy * (Ny - 2);

    % Determine z-slice range
    slice_z = params.slice_z;
    if slice_z > 0
        ss  = mk  * (slice_z - 1) + 1;
        se  = mk  * slice_z;
        ss2 = mk2 * (slice_z - 1) + 1;
        se2 = mk2 * slice_z;
    else
        ss  = 1;  se  = n_int;
        ss2 = 1;  se2 = n_int2;
    end

    % Compute B-field
    [~, ~, Bz_all] = evaluate_bfield(state_history, params);

    % Frame indices
    n_total = size(state_history, 2) - 1;
    frame_indices = floor(linspace(1, n_total, ip.Results.n_frames));

    fig = figure;
    fig.Position = [900, 200, 550, 800];
    fig.Color = [1, 1, 1];

    for idx = 1:length(frame_indices)
        t_idx = frame_indices(idx);

        % Extract psi for this time step
        psi_vec = state_history(1:n_int, t_idx);
        psi_slice = reshape(psi_vec(ss:se), Nx - 1, Ny - 1);

        % |ψ|² normalized
        psi_sq = abs(psi_slice).^2;
        psi_sq_norm = psi_sq ./ max(psi_sq(:), 1);

        % Bz slice
        bz_slice = reshape(real(Bz_all(ss2:se2, t_idx)), Nx - 2, Ny - 2);

        % Plot |ψ|²
        subplot(2, 1, 1)
        imagesc(x_coords, y_coords, psi_sq_norm)
        colorbar
        caxis([0, 1])
        axis equal tight off
        title('|\psi|^2')

        % Plot Bz
        subplot(2, 1, 2)
        imagesc(x_coords2, y_coords2, bz_slice)
        colorbar
        caxis(ip.Results.bz_clim)
        axis equal tight off
        title('B_z')

        sgtitle(sprintf('Applied Bz = %.4f, t = %.4f', ...
            params.B_real_z_history(t_idx), params.t(t_idx)))

        drawnow

        % Save GIF
        if ip.Results.save_gif
            if idx == 1
                gif_dir = ip.Results.gif_dir;
                if ~exist(gif_dir, 'dir')
                    mkdir(gif_dir);
                end
                filename = fullfile(gif_dir, ['trapGif', datestr(now, 30), '.gif']);
                gif(filename);
                save(strrep(filename, '.gif', '.mat'), 'params');
            else
                gif
            end
        end
    end
end
