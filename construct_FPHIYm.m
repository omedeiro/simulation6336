function FPHIY = construct_FPHIYm(x, y1, y2, y3, p)
    hx = p.hx;
%     hy = p.hy;
    hz = p.hz; 
    kappa = p.kappa; 
%     Nx = p.Nx;
%     Ny = p.Ny;
%     Nz = p.Nz;
    
    mk = (p.Nx+1)*(p.Ny+1);
    mj = (p.Nx+1);
    m = p.M2;
    M = p.m2;  

%     FPHIY = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

    FPHIY = sparse(M, 1, (kappa^2/hz^2).*(-y3(m+mj)+y3(m)+y3(m+mj-mk)-y3(m-mk))...
                    +(kappa^2/hx^2).*(-y1(m+mj)+y1(m)+y1(m-1+mj)-y1(m-1))...
                    +imag(exp(-1i*y2(m)).*conj(x(m)).*x(m+mj)) );

        
end                            

