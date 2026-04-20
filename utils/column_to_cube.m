function cube = column_to_cube(column, Nx, Ny, Nz)
% COLUMN_TO_CUBE  Reshape a column vector into a 3D array.
%
%   cube = column_to_cube(column, Nx, Ny, Nz)
%
%   Reshapes a column vector of length Nx*Ny*Nz into a 3D array of size
%   (Nx × Ny × Nz) using MATLAB's default column-major ordering.
%
%   Inputs:
%     column : double vector of length Nx*Ny*Nz
%     Nx     : number of points along x (first dimension)
%     Ny     : number of points along y (second dimension)
%     Nz     : number of points along z (third dimension)
%
%   See also CUBE_TO_COLUMN

    assert(numel(column) == Nx * Ny * Nz, ...
        'column_to_cube:sizeMismatch', ...
        'Length of column (%d) does not match Nx*Ny*Nz = %d.', ...
        numel(column), Nx * Ny * Nz);

    cube = reshape(column, [Nx, Ny, Nz]);
end
