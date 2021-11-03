% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)
% 
% x(i) = x(i+1)*exp(-1i*x(N+i)); %psi1
% x(N) = x(i-1)*exp(1i*x(N+i-1)); %psiN

% function F = analytical_f_xyz(x_int, y1_int, y2_int, y3_int, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
% function F = analytical_f_xyz(X, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
function F = eval_f(X, p, u)

    colN = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
    x_int = X(1:colN);
    y1_int = X(colN+1:2*colN);
    y2_int = X(2*colN+1:3*colN);
    y3_int = X(3*colN+1:4*colN);

    %%%%%%%% _int has size (p.Nx-1, p.Ny-1, p.Nz-1)
    x = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y1 = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y2 = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y3 = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    
    % change name for smeplicity
    M2 = p.M2;
    m2 = p.m2;
    m = p.m;
    
    m1x = p.m1x;
    m2x = p.m2x;
    mNx = p.mNx;
    mNxp1 = p.mNxp1;
    m1x_int = p.m1x_int;
    m2x_int = p.m2x_int;
    mNx_int = p.mNx_int;
    mNxp1_int = p.mNxp1_int;
    
    m1y = p.m1y;
    m2y = p.m2y;
    mNy = p.mNy;
    mNyp1 = p.mNyp1;
    m1y_int = p.m1y_int;
    m2y_int = p.m2y_int;
    mNy_int = p.mNy_int;
    mNyp1_int = p.mNyp1_int;
    
    m1z = p.m1z;
    m2z = p.m2z;
    mNz = p.mNz;
    mNzp1 = p.mNzp1;
    m1z_int = p.m1z_int;
    m2z_int = p.m2z_int;
    mNz_int = p.mNz_int;
    mNzp1_int = p.mNzp1_int;
    
    
    %%%%%%% set internal values
    x(M2) = x_int(m2);
    y1(M2) = y1_int(m2);
    y2(M2) = y2_int(m2);
    y3(M2) = y3_int(m2);


    %%% BOUNDARY CONDITIONS %%%
    
    mk = (p.Nx+1)*(p.Ny+1);
    mj = (p.Nx+1);
    
    % x boundaries
    periodic_x = 0;
            
            if periodic_x
                % periodic boundaries 
                x(m1x) = x(mNx); %34
                x(mNxp1) = x(m2x); %34

                y1(m1x) = y1(mNx); %34
                y1(Nxp1) = y1(m2x); %34
            else 
                % zero current on x
                x(m1x) = x(m2x).*exp(-1i*y1(m1x)); %35
                x(mNxp1) = x(mNx).*exp(1i*y1(mNx)); %35 
            end
           
            % Magnetic field x boundary conditions eq 37 
            
                y2(m1x_int) = -u.Bz*p.hx*p.hy + y1(m1x_int) - y1(m1x_int+mj) + y2(m2x_int);
                y2(mNxp1_int) = -u.Bz*p.hx*p.hy + y1(mNxp1_int) - y1(mNxp1_int+mj) + y2(mNx_int);

                y3(m1x_int) = u.By*p.hz*p.hx + y3(m2x_int) +y1(m1x_int) - y1(m1x_int+mk);
                y3(mNxp1_int) = u.By*p.hz*p.hx + y3(mNx_int) +y1(mNxp1_int) - y1(mNxp1_int+mk);
            
 
    % y boundaries
    periodic_y = 0;
            
            if periodic_y
                % periodic boundaries 
                x(m1y) = x(mNy); %34
                x(mNyp1) = x(m2y); %34

                y2(m1y) = y2(mNy); %34
                y2(Nyp1) = y2(m2y); %34
            else 
                % zero current on y
                x(m1y) = x(m2y).*exp(-1i*y2(m1y)); %35
                x(mNyp1) = x(mNy).*exp(1i*y2(mNy)); %35
            end

            % Magnetic field y boundary conditions eq 37 
            
                y1(m1y_int) = u.Bz*p.hx*p.hy + y1(m2y_int) + y2(m1y_int) - y2(m1y_int+1);
                y1(mNyp1_int) = u.Bz*p.hx*p.hy + y1(mNy_int) - y2(mNyp1_int) + y2(mNyp1_int+1);

                y3(m1y_int) = -u.Bx*p.hy*p.hz + y2(m1y_int) - y2(m1y_int+mk) + y3(m2y_int);
                y3(mNyp1_int) = -u.Bx*p.hy*p.hz + y2(mNyp1_int) - y2(mNyp1_int+mk) + y3(mNy_int);
            

    
     % z boundaries
    periodic_z = 0;

            if periodic_z
                % periodic boundaries 
                x(m1z) = x(mNz); %34
                x(mNzp1) = x(m2z); %34
                %%%%%%%%%%%%%%%%%%%%%%% WHY y1 periodic instead of y3, maybe error in paper?
                %%%%%%%%%%%%%%%%%%%%%%% y1(m1) is defined below
                y3(m1z) = y3(mNz); %34
                y3(Nzp1) = y3(m2z); %34
            else 
                % zero current on z
                x(m1z) = x(m2z).*exp(-1i*y3(m1z)); %35
                x(mNzp1) = x(mNz).*exp(1i*y3(mNz)); %35
            end
            
            % Magnetic field z boundary conditions eq 37 
            
                y1(m1z_int) = -u.By*p.hz*p.hx + y3(m1z_int) - y3(m1z_int+1) + y1(m2z_int);
                y1(mNzp1_int) = -u.By*p.hz*p.hx + y3(mNzp1_int) - y3(mNzp1_int+1) + y1(mNz_int);

                y2(m1z_int) = u.Bx*p.hy*p.hz + y2(m2z_int) + y3(m1z_int) - y3(m1z_int+mj);
                y2(mNzp1_int) = u.Bx*p.hy*p.hz + y2(mNz_int) + y3(mNzp1_int) - y3(mNzp1_int+mj);
            
            
  

    
    LPSIX = construct_LPSIXm(y1, p);
    LPSIY = construct_LPSIYm(y2, p);
    LPSIZ = construct_LPSIZm(y3, p);
      
    LPHIX = construct_LPHIXm(p);
    LPHIY = construct_LPHIYm(p);
    LPHIZ = construct_LPHIZm(p);
    
    FPSI = construct_FPSIm(x, p);
    FPHIX = construct_FPHIXm(x, y1, y2, y3, p);
    FPHIY = construct_FPHIYm(x, y1, y2, y3, p);
    FPHIZ = construct_FPHIZm(x, y1, y2, y3, p);
    
    % remove boundary rows (zeros) - NO EQUATIONS AT BOUNDARY
    
    LPSIX(~any(LPSIX,2), :) = [];
    LPSIY(~any(LPSIY,2), :) = [];
    LPSIZ(~any(LPSIZ,2), :) = [];

    LPHIX(~any(LPHIX,2), :) = [];
    LPHIY(~any(LPHIY,2), :) = [];
    LPHIZ(~any(LPHIZ,2), :) = [];

    %
    D = p.kappa^2;

    dPsidt = D*(LPSIX/p.hx^2 + LPSIY/p.hy^2 + LPSIZ/p.hz^2)*x + FPSI;
    dPhidtX = D*(LPHIY./p.hy^2 + LPHIZ./p.hz^2)*y1 + FPHIX;
    dPhidtY = D*(LPHIX./p.hx^2 + LPHIZ./p.hz^2)*y2 + FPHIY;
    dPhidtZ = D*(LPHIX./p.hx^2 + LPHIY./p.hy^2)*y3 + FPHIZ;

    F = [dPsidt; dPhidtX; dPhidtY; dPhidtZ];
end