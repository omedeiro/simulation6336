function F = eval_f(X, p, u)

    colN = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
    dim_x = (p.Nx+1)*(p.Ny+1)*(p.Nz+1);

    x_int = X(1:colN);
    y1_int = X(colN+1:2*colN);
    y2_int = X(2*colN+1:3*colN);
    y3_int = X(3*colN+1:4*colN);
    
    %%%%%%% set internal values
    x = sparse(p.M2, 1, x_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y1 = sparse(p.M2, 1, y1_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y2 = sparse(p.M2, 1, y2_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    y3 = sparse(p.M2, 1, y3_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    
    
    %% %%% BOUNDARY CONDITIONS %%%
        
    % x boundaries            
            if p.periodic_x
                % periodic boundaries 
                x = x + sparse(p.m1x_int, 1, x(p.mNx_int), dim_x, 1);
                x = x + sparse(p.mNxp1_int, 1, x(p.m2x_int), dim_x, 1);

                y1 = y1 + sparse(p.m1x_int, 1, y1(p.mNx_int), dim_x, 1);
                y1 = y1 + sparse(p.mNxp1_int, 1, y1(p.m2x_int), dim_x, 1);
            else 
                % zero current on x
                x = x + sparse(p.m1x_int, 1, x(p.m2x_int).*exp(-1i*y1(p.m1x_int)), dim_x, 1);
                x = x + sparse(p.mNxp1_int, 1, x(p.mNx_int).*exp(1i*y1(p.mNx_int)), dim_x, 1);
            end
            % Magnetic field x boundary conditions eq 37 
            y2 = y2 + sparse(p.m1x_int, 1, -u.Bz(p.m1x_int)*p.hx*p.hy + y2(p.m2x_int), dim_x, 1);
            y2 = y2 + sparse(p.mNxp1_int, 1, u.Bz(p.mNxp1_int)*p.hx*p.hy + y2(p.mNx_int), dim_x, 1);

            y3 = y3 + sparse(p.m1x_int, 1, u.By(p.m1x_int)*p.hz*p.hx + y3(p.m2x_int), dim_x, 1);
            y3 = y3 + sparse(p.mNxp1_int, 1, -u.By(p.mNxp1_int)*p.hz*p.hx + y3(p.mNx_int), dim_x, 1);
  
%     %%
    
%     [yy, xx, zz] = meshgrid(p.hx:p.hx:p.hx*(p.Nx+1), p.hy:p.hy:p.hy*(p.Ny+1), p.hz:p.hz:p.hz*(p.Nz+1));
% 
%     [yy2, xx2, zz2] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-2), p.hy:p.hy:p.hy*(p.Ny-2), p.hz:p.hz:p.hz*(p.Nz-2));
%     xx = cube2column(xx);
%     yy = cube2column(yy);
%     zz = cube2column(zz);
%     xx2 = cube2column(xx2);
%     yy2 = cube2column(yy2);
%     zz2 = cube2column(zz2);
%     figure
%     subplot(1,3,1)
%     scatter3(xx,yy,zz,36,abs(y1), 'filled', 'MarkerFaceAlpha', 0.5)
%     subplot(1,3,2)
%     scatter3(xx,yy,zz,36,abs(y2), 'filled', 'MarkerFaceAlpha', 0.5)
%     subplot(1,3,3)
%     scatter3(xx,yy,zz,36,abs(y3), 'filled', 'MarkerFaceAlpha', 0.5)
 
    % y boundaries            
            if p.periodic_y
                % periodic boundaries 
                x = x + sparse(p.m1y_int, 1, x(p.mNy_int), dim_x, 1);
                x = x + sparse(p.mNyp1_int, 1, x(p.m2y_int), dim_x, 1);

                y2 = y2 + sparse(p.m1y_int, 1, y2(p.mNy_int), dim_x, 1);
                y2 = y2 + sparse(p.mNyp1_int, 1, y2(p.m2y_int), dim_x, 1);
            else 
                % zero current on y
                x = x + sparse(p.m1y_int, 1, x(p.m2y_int).*exp(-1i*y2(p.m1y_int)), dim_x, 1);
                x = x + sparse(p.mNyp1_int, 1, x(p.mNy_int).*exp(1i*y2(p.mNy_int)), dim_x, 1);
            end

            % Magnetic field y boundary conditions eq 37 
            y1 = y1 + sparse(p.m1y_int, 1, u.Bz(p.m1y_int)*p.hx*p.hy + y1(p.m2y_int), dim_x, 1);
            y1 = y1 + sparse(p.mNyp1_int, 1, -u.Bz(p.mNyp1_int)*p.hx*p.hy + y1(p.mNy_int), dim_x, 1);
 
            y3 = y3 + sparse(p.m1y_int, 1, -u.Bx(p.m1y_int)*p.hy*p.hz + y3(p.m2y_int), dim_x, 1);
            y3 = y3 + sparse(p.mNyp1_int, 1, u.Bx(p.mNyp1_int)*p.hy*p.hz + y3(p.mNy_int), dim_x, 1);
                    
    % z boundaries
            if p.periodic_z
                % periodic boundaries 
                x = x + sparse(p.m1z_int, 1, x(p.mNz_int), dim_x, 1);
                x = x + sparse(p.mNzp1_int, 1, x(p.m2z_int), dim_x, 1);

                y3 = y3 + sparse(p.m1z_int, 1, y2(p.mNz_int), dim_x, 1);
                y3 = y3 + sparse(p.mNzp1_int, 1, y2(p.m2z_int), dim_x, 1);
            else 
                % zero current on z
                x = x + sparse(p.m1z_int, 1, x(p.m2z_int).*exp(-1i*y3(p.m1z_int)), dim_x, 1);
                x = x + sparse(p.mNzp1_int, 1, x(p.mNz_int).*exp(1i*y3(p.mNz_int)), dim_x, 1);
            end
            
            % Magnetic field z boundary conditions eq 37 
            y1 = y1 + sparse(p.m1z_int, 1, -u.By(p.m1z_int)*p.hz*p.hx + y1(p.m2z_int), dim_x, 1);
            y1 = y1 + sparse(p.mNzp1_int, 1, u.By(p.mNzp1_int)*p.hz*p.hx + y1(p.mNz_int), dim_x, 1);

            y2 = y2 + sparse(p.m1z_int, 1, u.Bx(p.m1z_int)*p.hy*p.hz + y2(p.m2z_int), dim_x, 1);
            y2 = y2 + sparse(p.mNzp1_int, 1, -u.Bx(p.mNzp1_int)*p.hy*p.hz + y2(p.mNz_int), dim_x, 1);

    


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