function F_PSI = construct_FPSIm(x, p)
    fun = @(a) (1-conj(a).*a).*a;
    
    F_PSI = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
    
    mk = (p.Nx+1)*(p.Ny+1);
    mj = (p.Nx+1);
    m = p.M2;
    M = p.m2;  

    F_PSI(M) = fun(x(m));


end
