% SETUP_PATHS  Add all simulation subdirectories to the MATLAB path.
%
%   Run this script from the project root or from any subdirectory.
%   It adds: core/, operators/, physics/, solvers/, visualization/,
%   analysis/, tests/, examples/, utils/, gui/.
%
%   Usage:
%     >> run('setup_paths.m')
%     -- or --
%     >> setup_paths
%
%   See also ADDPATH

    root = fileparts(mfilename('fullpath'));

    subdirs = {'core', 'operators', 'physics', 'solvers', ...
               'visualization', 'analysis', 'tests', 'examples', ...
               'utils', 'gui'};

    for i = 1:numel(subdirs)
        d = fullfile(root, subdirs{i});
        if exist(d, 'dir')
            addpath(d);
        end
    end

    v = tdgl_version();
    fprintf('TDGL 3D Simulator v%s — paths configured.\n', v);
