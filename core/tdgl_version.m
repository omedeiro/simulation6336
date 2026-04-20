function v = tdgl_version()
% TDGL_VERSION  Return the current version string of the TDGL package.
%
%   v = tdgl_version()
%
%   Reads the version from the VERSION file in the project root.
%   Falls back to a hardcoded value if the file is missing.
%
%   Example:
%     >> tdgl_version()
%     ans = '0.2.0'
%
%   See also SETUP_PATHS

    version_file = fullfile(fileparts(mfilename('fullpath')), '..', 'VERSION');
    if exist(version_file, 'file')
        fid = fopen(version_file, 'r');
        v = strtrim(fgetl(fid));
        fclose(fid);
    else
        v = '0.2.0';  % fallback
    end
end
