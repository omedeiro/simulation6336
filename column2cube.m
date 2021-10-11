function xcube = column2cube(X, Nx, Ny, Nz)
    xcube = zeros(Nx, Ny, Nz);
    i=1;
    j=1;
    k=1;
    for p = 1:numel(X)
        xcube(i,j,k) = X(p);
        if i==Nx && j~=Ny 
            i = 1;
            j = j+1;
        elseif i == Nx && j == Ny
            i = 1;
            j = 1;
            k = k + 1;
        else
            i = i + 1;
        end 
    end
end