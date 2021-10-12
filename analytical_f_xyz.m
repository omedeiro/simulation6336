% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)
% 
% x(i) = x(i+1)*exp(-1i*x(N+i)); %psi1
% x(N) = x(i-1)*exp(1i*x(N+i-1)); %psiN

% function F = analytical_f_xyz(x_int, y1_int, y2_int, y3_int, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
function F = analytical_f_xyz(X, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
    colN = (Nx-1)*(Ny-1)*(Nz-1);
    x_int = column2cube(X(1:colN), Nx-1, Ny-1, Nz-1);
    y1_int = column2cube(X(colN+1:2*colN), Nx-1, Ny-1, Nz-1);
    y2_int = column2cube(X(2*colN+1:3*colN), Nx-1, Ny-1, Nz-1);
    y3_int = column2cube(X(3*colN+1:4*colN), Nx-1, Ny-1, Nz-1);

    %%%%%%%% _int has size (Nx-1, Ny-1, Nz-1)
    x = zeros(Nx+1, Ny+1, Nz+1);
    y1 = zeros(Nx+1, Ny+1, Nz+1);
    y2 = zeros(Nx+1, Ny+1, Nz+1); 
    y3 = zeros(Nx+1, Ny+1, Nz+1);
    
    x(2:Nx, 2:Ny, 2:Nz) = x_int;
    y1(2:Nx, 2:Ny, 2:Nz) = y1_int;
    y2(2:Nx, 2:Ny, 2:Nz) = y2_int;
    y3(2:Nx, 2:Ny, 2:Nz) = y3_int;
    
    dPsidt = zeros(Nx+1, Ny+1, Nz+1);
    dPhidtX = zeros(Nx+1, Ny+1, Nz+1);
    dPhidtY = zeros(Nx+1, Ny+1, Nz+1);
    dPhidtZ = zeros(Nx+1, Ny+1, Nz+1);
    
    D = kappa^2;

    %%% BOUNDARY CONDITIONS %%%
    
    x(:,:,1) = x(:,:,Nz); %34
    x(:,:,Nz+1) = x(:,:,2); %34
    
    y1(:,:,1) = y1(:,:,Nz); %34
    y1(:,:,Nz+1) = y1(:,:,2); %34
    
    x(1,:,:) = x(2,:,:).*exp(-1i*y1(1,:,:)); %35
    x(Nx+1,:,:) = x(Nx,:,:).*exp(1i*y1(Nx,:,:)); %35    
    
    x(:,1,:) = x(:,2,:).*exp(-1i*y2(:,1,:)); %35
    x(:,Ny+1,:) = x(:,Ny,:).*exp(1i*y2(:,Ny,:)); %35
    
    % BCs on yz for Bx (36) 
    om=1;
    if om ~= 1
        if Bx==0
            row_B = ones(1, Ny+1);
        else
    %         row_B = [1 : Bx*hy*hz : Bx*hy*hz*(Ny+1)];
            row_B = linspace(1, Bx*hy*hz*(Ny+1), Ny+1);
        end
        y3_b = row_B'*ones(1,Nz+1);

        %%%% phi_z in first layer yz
        % first face
        y3(1, :, :) = y3_b;
        y2(1, :, :) = ones(Ny+1, Nz+1);
    %     % second face
    %     y3(2, :, :) = y3_b;
    %     y2(2, :, :) = zeros(Ny+1, Nz+1);
        %connections along x
        y1(1, :, :) = ones(Ny+1, Nz+1);

        %%%% phi_z in last layer yz
        % first face
        y3(Nx+1, :, :) = y3_b;
        y2(Nx+1, :, :) = ones(Ny+1, Nz+1);
    %     % second face
    %     y3(Nx+1, :, :) = y3_b;
    %     y2(Nx+1, :, :) = zeros(Ny+1, Nz+1);

        %connections along x
        y1(Nx+1, :, :) = ones(Ny+1, Nz+1);

    else
        % Magnetic field boundary conditions eq 37
        By=0; % we can update later
        Bz=0; % we can update later
        for kk = 2:Nz
            for ii = 2:Nx
                y1(ii,1,kk) = -Bz*hx*hy + y1(ii,2,kk) + y2(ii,1,kk) - y2(ii+1,1,kk);
                y1(ii,Ny+1,kk) = Bz*hx*hy + y1(ii,Ny,kk) - y2(ii,Ny+1,kk) + y2(ii+1,Ny+1,kk);
                
                y2(ii,1,kk) = -Bz*hx*hy - y2(ii+1,1,kk) + y1(ii,1,kk) - y1(ii,2,kk);
                y2(ii,Ny+1,kk) = -Bz*hx*hy - y2(ii+1,Ny+1,kk) + y1(ii,Ny+1,kk) - y1(ii,Ny,kk);
                
                y3(ii,1,kk) = -Bx*hy*hz + y2(ii,1,kk) - y2(ii,1,kk+1) + y3(ii,2,kk);
                y3(ii,Ny+1,kk) = Bx*hy*hz - y2(ii,Ny+1,kk) + y2(ii,Ny+1,kk+1) + y3(ii,Ny,kk);
            end
        end
        for jj = 2:Ny
            for ii = 2:Nx
                y2(ii,jj,1) = -Bx*hy*hz + y2(ii,jj,2) + y3(ii,jj,1) - y3(ii,jj+1,1);
                y2(ii,jj,Nz+1) = Bx*hy*hz + y2(ii,jj,Nz) - y3(ii,jj,Nz+1) + y3(ii,jj+1,Nz+1);

                y3(ii,jj,1) = -Bx*hy*hz - y3(ii,jj+1,1) + y2(ii,jj,1) - y2(ii,jj,2);
                y3(ii,jj,Nz+1) = -Bx*hy*hz - y3(ii,jj+1,Nz) + y2(ii,jj,Nz+1) - y2(ii,jj,Nz);           
            end
        end
        for kk = 2:Nz
            for jj = 2:Ny
                y3(1,jj,kk) = -By*hx*hz + y3(2,jj,kk) + y1(1,jj,kk) - y1(1,jj,kk+1);
                y3(Nx+1,jj,kk) = By*hy*hz + y3(Nx,jj,kk) - y1(Nx+1,jj,kk) + y1(Nx+1,jj,kk+1);

                y1(1,jj,kk) = -By*hx*hz - y1(1,jj,kk+1) + y3(1,jj,kk) - y3(2,jj,kk);
                y1(Nx+1,jj,kk) = -By*hy*hz - y1(Nx,jj,kk+1) + y3(Nx+1,jj,kk) - y3(Nx,jj,kk);
                
                y2(1,jj,kk) = -Bz*hy*hx + y1(1,jj,kk) - y1(1,jj+1,kk) + y2(2,jj,kk);
                y2(Nx+1,jj,kk) = Bz*hy*hx - y1(Nx+1,jj,kk) + y1(Nx+1,jj+1,kk) + y2(Nx,jj,kk);
            end
        end

    end
    
    LPSIX = construct_LPSIX(y1, Nx, Ny, Nz);
    LPSIY = construct_LPSIY(y2, Nx, Ny, Nz);
    LPSIZ = construct_LPSIZ(y3, Nx, Ny, Nz);
      
    LPHIX = construct_LPHIX(hx, hy, hz, kappa, Nx, Ny, Nz);
    LPHIY = construct_LPHIX(hx, hy, hz, kappa, Nx, Ny, Nz);
    LPHIZ = construct_LPHIX(hx, hy, hz, kappa, Nx, Ny, Nz);
    
    FPSI = construct_FPSI(x, Nx, Ny, Nz);
    FPHIX = construct_FPHIX(x, y1, y2, y3, hy, hz, kappa, Nx, Ny, Nz);
    FPHIY = construct_FPHIY(x, y1, y2, y3, hy, hz, kappa, Nx, Ny, Nz);
    FPHIZ = construct_FPHIZ(x, y1, y2, y3, hy, hz, kappa, Nx, Ny, Nz);
    
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
    dPsidt = D*(LPSIX/hx^2 + LPSIY/hy^2 + LPSIZ/hz^2)*u_x + FPSI;
    dPhidtX = D*(0/hx^2 + LPHIY/hy^2 + LPHIZ/hz^2)*u_y1 + FPHIX;
    dPhidtY = D*(LPHIX/hx^2 + 0/hy^2 + LPHIZ/hz^2)*u_y2 + FPHIY;
    dPhidtZ = D*(LPHIX/hx^2 + LPHIY/hy^2 + 0/hz^2)*u_y3 + FPHIZ;
    


    F = [dPsidt; dPhidtX; dPhidtY; dPhidtZ];
end