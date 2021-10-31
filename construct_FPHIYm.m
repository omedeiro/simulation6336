function FPHIY = construct_FPHIYm(X, p)
    hx = p.hx;
    hy = p.hy;
    hz = p.hz; 
    kappa = p.kappa; 
    Nx = p.Nx;
    Ny = p.Ny;
    Nz = p.Nz;
    
    colN = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);
    x = X(1:colN);
    y1 = X(colN+1:2*colN);
    y2 = X(2*colN+1:3*colN);
    y3 = X(3*colN+1:4*colN);

    FPHIY = sparse((Nx-1)*(Ny-1)*(Nz-1),1);
    for k = 2 : Nz
        for j = 2 : Ny
            for i = 2 : Nx
                mj = (Nx+1)*j-(Nx+1)*(j-1);
                mk = (Nx+1)*(Ny+1)*k-(Nx+1)*(Ny+1)*(k-1);
                m = i + (Nx+1)*(j-1)+(Nx+1)*(Ny+1)*(k-1);
                M = i + (Nx-1)*(j-1)+(Nx-1)*(Ny-1)*(k-1);
                
                FPHIY(M) = (kappa^2/hz^2)*(-y3(m+mj)+y3(m)+y3(m+mj-mk)-y3(m-mk))...
                                +(kappa^2/hx^2)*(-y1(m+mj)+y1(m)+y1(m-1+mj)-y1(m-1))...
                                +imag(exp(-1i*y2(m))*conj(x(m))*x(m+mj));
            
            end
        end
    end
end                            

