function F = eval_f(X, p, u)

    colN = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
    x_int = X(1:colN);
    y1_int = X(colN+1:2*colN);
    y2_int = X(2*colN+1:3*colN);
    y3_int = X(3*colN+1:4*colN);

    % change name for simplicity 
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
    x = sparse(M2, 1, x_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y1 = sparse(M2, 1, y1_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y2 = sparse(M2, 1, y2_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y3 = sparse(M2, 1, y3_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    
    %%% BOUNDARY CONDITIONS %%%
    
    dim_x = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);
    
    % x boundaries
    periodic_x = p.periodic_x;
            
            if periodic_x
                % periodic boundaries 
                x = x + sparse(m1x_int, 1, x(mNx_int), dim_x, 1);
                x = x + sparse(mNxp1_int, 1, x(m2x_int), dim_x, 1);

                y1 = y1 + sparse(m1x_int, 1, y1(mNx_int), dim_x, 1);
                y1 = y1 + sparse(mNxp1_int, 1, y1(m2x_int), dim_x, 1);
            else 
                % zero current on x
                x = x + sparse(m1x_int, 1, x(m2x_int).*exp(-1i*y1(m1x_int)), dim_x, 1);
                x = x + sparse(mNxp1_int, 1, x(mNx_int).*exp(1i*y1(mNx_int)), dim_x, 1);
            end
            % Magnetic field x boundary conditions eq 37 
            y2 = y2 + sparse(m1x_int, 1, -u.Bz(m1x_int)*p.hx*p.hy + y2(m2x_int), dim_x, 1);
            y2 = y2 + sparse(mNxp1_int, 1, u.Bz(mNxp1_int)*p.hx*p.hy + y2(mNx_int), dim_x, 1);

            y3 = y3 + sparse(m1x_int, 1, u.By(m1x_int)*p.hz*p.hx + y3(m2x_int), dim_x, 1);
            y3 = y3 + sparse(mNxp1_int, 1, -u.By(mNxp1_int)*p.hz*p.hx + y3(mNx_int), dim_x, 1);
                
 
    % y boundaries
    periodic_y = p.periodic_y;
            
            if periodic_y
                % periodic boundaries 
                x = x + sparse(m1y_int, 1, x(mNy_int), dim_x, 1);
                x = x + sparse(mNyp1_int, 1, x(m2y_int), dim_x, 1);

                y2 = y2 + sparse(m1y_int, 1, y2(mNy_int), dim_x, 1);
                y2 = y2 + sparse(mNyp1_int, 1, y2(m2y_int), dim_x, 1);
            else 
                % zero current on y
                x = x + sparse(m1y_int, 1, x(m2y_int).*exp(-1i*y2(m1y_int)), dim_x, 1);
                x = x + sparse(mNyp1_int, 1, x(mNy_int).*exp(1i*y2(mNy_int)), dim_x, 1);
            end

            % Magnetic field y boundary conditions eq 37 
            y1 = y1 + sparse(m1y_int, 1, u.Bz(m1y_int)*p.hx*p.hy + y1(m2y_int), dim_x, 1);
            y1 = y1 + sparse(mNyp1_int, 1, -u.Bz(mNyp1_int)*p.hx*p.hy + y1(mNy_int), dim_x, 1);
 
            y3 = y3 + sparse(m1y_int, 1, -u.Bx(m1y_int)*p.hy*p.hz + y3(m2y_int), dim_x, 1);
            y3 = y3 + sparse(mNyp1_int, 1, u.Bx(mNyp1_int)*p.hy*p.hz + y3(mNy_int), dim_x, 1);
                    
    % z boundaries
    periodic_z = p.periodic_z;

            if periodic_z
                % periodic boundaries 
                x = x + sparse(m1z_int, 1, x(mNz_int), dim_x, 1);
                x = x + sparse(mNzp1_int, 1, x(m2z_int), dim_x, 1);

                y3 = y3 + sparse(m1z_int, 1, y2(mNz_int), dim_x, 1);
                y3 = y3 + sparse(mNzp1_int, 1, y2(m2z_int), dim_x, 1);
            else 
                % zero current on z
                x = x + sparse(m1z_int, 1, x(m2z_int).*exp(-1i*y3(m1z_int)), dim_x, 1);
                x = x + sparse(mNzp1_int, 1, x(mNz_int).*exp(1i*y3(mNz_int)), dim_x, 1);
            end
            
            % Magnetic field z boundary conditions eq 37 
            y1 = y1 + sparse(m1z_int, 1, -u.By(m1z_int)*p.hz*p.hx + y1(m2z_int), dim_x, 1);
            y1 = y1 + sparse(mNzp1_int, 1, u.By(mNzp1_int)*p.hz*p.hx + y1(mNz_int), dim_x, 1);

            y2 = y2 + sparse(m1z_int, 1, u.Bx(m1z_int)*p.hy*p.hz + y2(m2z_int), dim_x, 1);
            y2 = y2 + sparse(mNzp1_int, 1, -u.Bx(mNzp1_int)*p.hy*p.hz + y2(mNz_int), dim_x, 1);

    
%     %%
%     
%     [yy, xx, zz] = meshgrid(p.hx:p.hx:p.hx*(p.Nx+1), p.hy:p.hy:p.hy*(p.Ny+1), p.hz:p.hz:p.hz*(p.Nz+1));
% 
%     [yy2, xx2, zz2] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-2), p.hy:p.hy:p.hy*(p.Ny-2), p.hz:p.hz:p.hz*(p.Nz-2));
%     xx = cube2column(xx);
%     yy = cube2column(yy);
%     zz = cube2column(zz);
%     xx2 = cube2column(xx2);
%     yy2 = cube2column(yy2);
%     zz2 = cube2column(zz2);
%     figure(1)
%     subplot(1,3,1)
%     scatter3(xx,yy,zz,36,y1, 'filled', 'MarkerFaceAlpha', 0.5)
%     subplot(1,3,2)
%     scatter3(xx,yy,zz,36,y2, 'filled', 'MarkerFaceAlpha', 0.5)
%     subplot(1,3,3)
%     scatter3(xx,yy,zz,36,y3, 'filled', 'MarkerFaceAlpha', 0.5)

%%
    LPSIX = construct_LPSIXm(y1, p);
    LPSIY = construct_LPSIYm(y2, p);
    LPSIZ = construct_LPSIZm(y3, p);
    
    LPHIX = p.LPHIX;
    LPHIY = p.LPHIY;
    LPHIZ = p.LPHIZ;
    
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

    
    D = p.kappa^2;

    dPsidt = D*(LPSIX/p.hx^2 + LPSIY/p.hy^2 + LPSIZ/p.hz^2)*x + FPSI;
    dPhidtX = D*(LPHIY./p.hy^2 + LPHIZ./p.hz^2)*y1 + FPHIX;
    dPhidtY = D*(LPHIX./p.hx^2 + LPHIZ./p.hz^2)*y2 + FPHIY;
    dPhidtZ = D*(LPHIX./p.hx^2 + LPHIY./p.hy^2)*y3 + FPHIZ;

    F = [dPsidt; dPhidtX; dPhidtY; dPhidtZ];
end