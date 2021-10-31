function F_PSI = construct_FPSIm(x, p)
    fun = @(a) (1-conj(a).*a).*a;
    
    F_PSI = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
    
    for k = 2 : p.Nz
        for j = 2 : p.Ny
            for i = 2 : p.Nx
                m = i+(p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                M = i-1 + (p.Nx-1)*(j-2)+(p.Nx-1)*(p.Ny-1)*(k-2);
                
                F_PSI(M) = fun(x(m));
            end
        end
    end

end
