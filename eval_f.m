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
    
    %%%%%%% set internal values
    for k = 2 : p.Nz
        for j = 2 : p.Ny
            for i = 2 : p.Nx
                M = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                m = i-1 + (p.Nx-1)*(j-2)+(p.Nx-1)*(p.Ny-1)*(k-2);
    
                x(M) = x_int(m);
                y1(M) = y1_int(m);
                y2(M) = y2_int(m);
                y3(M) = y3_int(m);
            end
        end
    end

    %%% BOUNDARY CONDITIONS %%%
    
    mk = (p.Nx+1)*(p.Ny+1);
    mj = (p.Nx+1);
    
    % x boundaries
    periodic_x = 0;
    for k = 1 : p.Nz+1
        for j = 1 : p.Ny+1
            m1 = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            m2 = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            mNx = Nx + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            mNxp1 = Nx+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            
            if periodic_x
                % periodic boundaries 
                x(m1) = x(mNx); %34
                x(mNxp1) = x(m2); %34

                y1(m1) = y1(mNx); %34
                y1(Nxp1) = y1(m2); %34
            else 
                % zero current on x
                x(m1) = x(m2).*exp(-1i*y1(m1)); %35
                x(mNxp1) = x(mNx).*exp(1i*y1(mNx)); %35 
            end
           
            % Magnetic field x boundary conditions eq 37 
            if k ~= 1 && k ~= p.Nz+1 && j ~= 1 && j ~= p.Ny+1
                y2(m1) = -u.Bz*p.hx*p.hy + y1(m1) - y1(m1+mj) + y2(m2);
                y2(mNxp1) = -u.Bz*p.hx*p.hy + y1(mNxp1) - y1(mNxp1+mj) + y2(mNx);

                y3(m1) = u.By*p.hz*p.hx + y3(m2) +y1(m1) - y1(m1+mk);
                y3(mNxp1) = u.By*p.hz*p.hx + y3(mNx) +y1(mNxp1) - y1(mNxp1+mk);
            end
        end
    end
    
    % y boundaries
    periodic_y = 0;
    for k = 1 : p.Nz+1
        for i = 1 : p.Nx+1 
            m1 = i + (p.Nx+1)*(p.Ny+1)*(k-1);
            m2 = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            mNy = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            mNyp1 = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
            
            if periodic_y
                % periodic boundaries 
                x(m1) = x(mNy); %34
                x(mNyp1) = x(m2); %34

                y2(m1) = y2(mNy); %34
                y2(Nyp1) = y2(m2); %34
            else 
                % zero current on y
                x(m1) = x(m2).*exp(-1i*y2(m1)); %35
                x(mNyp1) = x(mNy).*exp(1i*y2(mNy)); %35
            end

            % Magnetic field y boundary conditions eq 37 
            if k ~= 1 && k ~= p.Nz+1 && i ~= 1 && i ~= p.Nx+1
                y1(m1) = u.Bz*p.hx*p.hy + y1(m2) + y2(m1) - y2(m1+1);
                y1(mNyp1) = u.Bz*p.hx*p.hy + y1(mNy) - y2(mNyp1) + y2(mNyp1+1);

                y3(m1) = -u.Bx*p.hy*p.hz + y2(m1) - y2(m1+mk) + y3(m2);
                y3(mNyp1) = -u.Bx*p.hy*p.hz + y2(mNyp1) - y2(mNyp1+mk) + y3(mNy);
            end
        end
    end
    
     % z boundaries
    periodic_z = 0;
    for j = 1 : p.Ny+1
        for i = 1 : p.Nx+1
            m1 = i + (p.Nx+1)*(j-1);
            m2 = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
            mNz = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
            mNzp1 = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*p.Nz;

            if periodic_z
                % periodic boundaries 
                x(m1) = x(mNz); %34
                x(mNzp1) = x(m2); %34
                %%%%%%%%%%%%%%%%%%%%%%% WHY y1 periodic instead of y3, maybe error in paper?
                %%%%%%%%%%%%%%%%%%%%%%% y1(m1) is defined below
                y3(m1) = y3(mNz); %34
                y3(Nzp1) = y3(m2); %34
            else 
                % zero current on z
                x(m1) = x(m2).*exp(-1i*y3(m1)); %35
                x(mNzp1) = x(mNz).*exp(1i*y3(mNz)); %35
            end
            
            % Magnetic field z boundary conditions eq 37 
            if j ~= 1 && j ~= p.Ny+1 && i ~= 1 && i ~= p.Nx+1
                y1(m1) = -u.By*p.hz*p.hx + y3(m1) - y3(m1+1) + y1(m2);
                y1(mNzp1) = -u.By*p.hz*p.hx + y3(mNzp1) - y3(mNzp1+1) + y1(mNz);

                y2(m1) = u.Bx*p.hy*p.hz + y2(m2) + y3(m1) - y3(m1+mj);
                y2(mNzp1) = u.Bx*p.hy*p.hz + y2(mNz) + y3(mNzp1) - y3(mNzp1+mj);
            end
            
        end
    end

    
    LPSIX = construct_LPSIX(y1, p.Nx, p.Ny, p.Nz);
    LPSIY = construct_LPSIY(y2, p.Nx, p.Ny, p.Nz);
    LPSIZ = construct_LPSIZ(y3, p.Nx, p.Ny, p.Nz);
      
    LPHIX = construct_LPHIX(p.hx, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    LPHIY = construct_LPHIY(p.hx, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    LPHIZ = construct_LPHIZ(p.hx, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    
    FPSI = construct_FPSI(x, p.Nx, p.Ny, p.Nz);
    FPHIX = construct_FPHIX(x, y1, y2, y3, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    FPHIY = construct_FPHIY(x, y1, y2, y3, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    FPHIZ = construct_FPHIZ(x, y1, y2, y3, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    
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