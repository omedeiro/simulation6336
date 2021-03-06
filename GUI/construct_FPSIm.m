function F_PSI = construct_FPSIm(x, p)
    fun = @(a) (1-conj(a).*a).*a;
    m = p.M2;
    M = p.m2; 
    
    %F_PSI = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
    
    F_PSI = sparse(M,1,fun(x(m)));
    


end
