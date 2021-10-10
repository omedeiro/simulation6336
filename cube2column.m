function ucol = cube2column(u)
    Nx = size(u, 1);
    Ny = size(u, 2);
    Nz = size(u, 3);
    ucol = zeros(Nx*Ny*Nz,1);
    for k = 1:Nz
        for j = 1:Ny
            for i = 1:Nx
                m = Nx*(j-1)+Nx*Ny*(k-1);
                ucol(i+m) = u(i,j,k);
            end
        end
end