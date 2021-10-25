% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)
% 
% x(i) = x(i+1)*exp(-1i*x(N+i)); %psi1
% x(N) = x(i-1)*exp(1i*x(N+i-1)); %psiN

% function F = analytical_f_xyz(x_int, y1_int, y2_int, y3_int, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
% function F = analytical_f_xyz(X, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
function F = analytical_f_xyz(X, p, u)

    colN = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
    x_int = column2cube(X(1:colN), p.Nx-1, p.Ny-1, p.Nz-1);
    y1_int = column2cube(X(colN+1:2*colN), p.Nx-1, p.Ny-1, p.Nz-1);
    y2_int = column2cube(X(2*colN+1:3*colN), p.Nx-1, p.Ny-1, p.Nz-1);
    y3_int = column2cube(X(3*colN+1:4*colN), p.Nx-1, p.Ny-1, p.Nz-1);

    %%%%%%%% _int has size (p.Nx-1, p.Ny-1, p.Nz-1)
    x = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    y1 = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    y2 = zeros(p.Nx+1, p.Ny+1, p.Nz+1); 
    y3 = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    
    x(2:p.Nx, 2:p.Ny, 2:p.Nz) = x_int;
    y1(2:p.Nx, 2:p.Ny, 2:p.Nz) = y1_int;
    y2(2:p.Nx, 2:p.Ny, 2:p.Nz) = y2_int;
    y3(2:p.Nx, 2:p.Ny, 2:p.Nz) = y3_int;
    
    dPsidt = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    dPhidtX = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    dPhidtY = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    dPhidtZ = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    
    D = p.kappa^2;

    %%% BOUNDARY CONDITIONS %%%
    
    x(:,:,1) = x(:,:,p.Nz); %34
    x(:,:,p.Nz+1) = x(:,:,2); %34
    
    y1(:,:,1) = y1(:,:,p.Nz); %34
    y1(:,:,p.Nz+1) = y1(:,:,2); %34
    
    x(1,:,:) = x(2,:,:).*exp(-1i*y1(1,:,:)); %35
    x(p.Nx+1,:,:) = x(p.Nx,:,:).*exp(1i*y1(p.Nx,:,:)); %35    
    
    x(:,1,:) = x(:,2,:).*exp(-1i*y2(:,1,:)); %35
    x(:,p.Ny+1,:) = x(:,p.Ny,:).*exp(1i*y2(:,p.Ny,:)); %35
    
    % BCs on yz for u.Bx (36) 
    om=1;
    if om ~= 1
        if u.Bx==0
            row_B = ones(1, p.Ny+1);
        else
    %         row_B = [1 : u.Bx*p.hy*p.hz : u.Bx*p.hy*p.hz*(p.Ny+1)];
            row_B = linspace(1, u.Bx*p.hy*p.hz*(p.Ny+1), p.Ny+1);
        end
        y3_b = row_B'*ones(1,p.Nz+1);

        %%%% phi_z in first layer yz
        % first face
        y3(1, :, :) = y3_b;
        y2(1, :, :) = ones(p.Ny+1, p.Nz+1);
    %     % second face
    %     y3(2, :, :) = y3_b;
    %     y2(2, :, :) = zeros(p.Ny+1, p.Nz+1);
        %connections along x
        y1(1, :, :) = ones(p.Ny+1, p.Nz+1);

        %%%% phi_z in last layer yz
        % first face
        y3(p.Nx+1, :, :) = y3_b;
        y2(p.Nx+1, :, :) = ones(p.Ny+1, p.Nz+1);
    %     % second face
    %     y3(p.Nx+1, :, :) = y3_b;
    %     y2(p.Nx+1, :, :) = zeros(p.Ny+1, p.Nz+1);

        %connections along x
        y1(p.Nx+1, :, :) = ones(p.Ny+1, p.Nz+1);

    else
        % Magnetic field boundary conditions eq 37
        u.By=0; % we can update later
        u.Bz=0; % we can update later
        for kk = 2:p.Nz
            for ii = 2:p.Nx
                y1(ii,1,kk) = -u.Bz*p.hx*p.hy + y1(ii,2,kk) + y2(ii,1,kk) - y2(ii+1,1,kk);
                y1(ii,p.Ny+1,kk) = u.Bz*p.hx*p.hy + y1(ii,p.Ny,kk) - y2(ii,p.Ny+1,kk) + y2(ii+1,p.Ny+1,kk);
                
                y2(ii,1,kk) = -u.Bz*p.hx*p.hy - y2(ii+1,1,kk) + y1(ii,1,kk) - y1(ii,2,kk);
                y2(ii,p.Ny+1,kk) = -u.Bz*p.hx*p.hy - y2(ii+1,p.Ny+1,kk) + y1(ii,p.Ny+1,kk) - y1(ii,p.Ny,kk);
                
                y3(ii,1,kk) = -u.Bx*p.hy*p.hz + y2(ii,1,kk) - y2(ii,1,kk+1) + y3(ii,2,kk);
                y3(ii,p.Ny+1,kk) = u.Bx*p.hy*p.hz - y2(ii,p.Ny+1,kk) + y2(ii,p.Ny+1,kk+1) + y3(ii,p.Ny,kk);
            end
        end
        for jj = 2:p.Ny
            for ii = 2:p.Nx
                y2(ii,jj,1) = -u.Bx*p.hy*p.hz + y2(ii,jj,2) + y3(ii,jj,1) - y3(ii,jj+1,1);
                y2(ii,jj,p.Nz+1) = u.Bx*p.hy*p.hz + y2(ii,jj,p.Nz) - y3(ii,jj,p.Nz+1) + y3(ii,jj+1,p.Nz+1);

                y3(ii,jj,1) = -u.Bx*p.hy*p.hz - y3(ii,jj+1,1) + y2(ii,jj,1) - y2(ii,jj,2);
                y3(ii,jj,p.Nz+1) = -u.Bx*p.hy*p.hz - y3(ii,jj+1,p.Nz) + y2(ii,jj,p.Nz+1) - y2(ii,jj,p.Nz);           
            end
        end
        for kk = 2:p.Nz
            for jj = 2:p.Ny
                y3(1,jj,kk) = -u.By*p.hx*p.hz + y3(2,jj,kk) + y1(1,jj,kk) - y1(1,jj,kk+1);
                y3(p.Nx+1,jj,kk) = u.By*p.hy*p.hz + y3(p.Nx,jj,kk) - y1(p.Nx+1,jj,kk) + y1(p.Nx+1,jj,kk+1);

                y1(1,jj,kk) = -u.By*p.hx*p.hz - y1(1,jj,kk+1) + y3(1,jj,kk) - y3(2,jj,kk);
                y1(p.Nx+1,jj,kk) = -u.By*p.hy*p.hz - y1(p.Nx,jj,kk+1) + y3(p.Nx+1,jj,kk) - y3(p.Nx,jj,kk);
                
                y2(1,jj,kk) = -u.Bz*p.hy*p.hx + y1(1,jj,kk) - y1(1,jj+1,kk) + y2(2,jj,kk);
                y2(p.Nx+1,jj,kk) = u.Bz*p.hy*p.hx - y1(p.Nx+1,jj,kk) + y1(p.Nx+1,jj+1,kk) + y2(p.Nx,jj,kk);
            end
        end

    end
    
    LPSIX = construct_LPSIX(y1, p.Nx, p.Ny, p.Nz);
    LPSIY = construct_LPSIY(y2, p.Nx, p.Ny, p.Nz);
    LPSIZ = construct_LPSIZ(y3, p.Nx, p.Ny, p.Nz);
      
    LPHIX = construct_LPHIX(p.hx, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    LPHIY = construct_LPHIX(p.hx, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    LPHIZ = construct_LPHIX(p.hx, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    
    FPSI = construct_FPSI(x, p.Nx, p.Ny, p.Nz);
    FPHIX = construct_FPHIX(x, y1, y2, y3, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    FPHIY = construct_FPHIY(x, y1, y2, y3, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    FPHIZ = construct_FPHIZ(x, y1, y2, y3, p.hy, p.hz, p.kappa, p.Nx, p.Ny, p.Nz);
    
    % create column vectors
    u_x = cube2column(x);
    u_y1 = cube2column(y1);
    u_y2 = cube2column(y2);
    u_y3 = cube2column(y3);
    
    % remove boundary rows (zeros) - NO EQUATIONS AT BOUNDARY
    
    LPSIX(~any(LPSIX,2), :) = [];
    LPSIY(~any(LPSIY,2), :) = [];
    LPSIZ(~any(LPSIZ,2), :) = [];

    LPHIX(~any(LPHIX,2), :) = [];
    LPHIY(~any(LPHIY,2), :) = [];
    LPHIZ(~any(LPHIZ,2), :) = [];

%     FPSI(~any(FPSI,2), :) = [];
%     FPHIX(~any(FPHIX,2), :) = [];
%     FPHIY(~any(FPHIY,2), :) = [];
%     FPHIZ(~any(FPHIZ,2), :) = [];
    
    %
    dPsidt = D*(LPSIX/p.hx^2 + LPSIY/p.hy^2 + LPSIZ/p.hz^2)*u_x + FPSI;
    dPhidtX = D*(0/p.hx^2 + LPHIY/p.hy^2 + LPHIZ/p.hz^2)*u_y1 + FPHIX;
    dPhidtY = D*(LPHIX/p.hx^2 + 0/p.hy^2 + LPHIZ/p.hz^2)*u_y2 + FPHIY;
    dPhidtZ = D*(LPHIX/p.hx^2 + LPHIY/p.hy^2 + 0/p.hz^2)*u_y3 + FPHIZ;
    


    F = [dPsidt; dPhidtX; dPhidtY; dPhidtZ];
%     F = abs(F);
end