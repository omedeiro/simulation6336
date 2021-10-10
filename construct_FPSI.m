function F_PSI = construct_FPSI(x, Nx, Ny, Nz)
    fun = @(a) (1-conj(a).*a).*a;
    F_PSI = cube2column(fun(  x(2:Nx,2:Ny,2:Nz)  ));
end
