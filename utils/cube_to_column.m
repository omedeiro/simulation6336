function column = cube_to_column(cube)
% CUBE_TO_COLUMN  Reshape a 3D array into a column vector.
%
%   column = cube_to_column(cube)
%
%   Flattens a 3D array (Nx × Ny × Nz) into a column vector of length
%   Nx*Ny*Nz using MATLAB's default column-major ordering:
%       column(m) = cube(i, j, k)
%   with m = i + (j-1)*Nx + (k-1)*Nx*Ny.
%
%   See also COLUMN_TO_CUBE

    column = cube(:);
end
