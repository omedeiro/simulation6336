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
    if Bx==0
        row_B = zeros(1, Ny+1);
    else
        row_B = [1 : Bx*hy*hz : Bx*hy*hz*(Ny+1)];
    end
    y3_b = row_B'*ones(1,Nz+1);

    %%%% phi_z in first layer yz
    % first face
    y3(1, :, :) = y3_b;
    y2(1, :, :) = zeros(Ny+1, Nz+1);
    % second face
    y3(2, :, :) = y3_b;
    y2(2, :, :) = zeros(Ny+1, Nz+1);
    %connections along x
    y1(1, :, :) = zeros(Ny+1, Nz+1);

    %%%% phi_z in last layer yz
    % first face
    y3(Nx, :, :) = y3_b;
    y2(Nx, :, :) = zeros(Ny+1, Nz+1);
    % second face
    y3(Nx+1, :, :) = y3_b;
    y2(Nx+1, :, :) = ones(Ny+1, Nz+1);
    %connections along x
    y1(Nx, :, :) = zeros(Ny+1, Nz+1);
    
    LPSIX = construct_LPSIX(y1, Nx, Ny, Nz);
    LPSIY = construct_LPSIY(y2, Nx, Ny, Nz);
    LPSIZ = construct_LPSIZ(y3, Nx, Ny, Nz);
      
    LPHIX = construct_LPHIX(hx, hy, hz, kappa, Nx, Ny, Nz);
    LPHIY = construct_LPHIX(hx, hy, hz, kappa, Nx, Ny, Nz);
    LPHIZ = construct_LPHIX(hx, hy, hz, kappa, Nx, Ny, Nz);
    
    FPSI = construct_FPSI(x);
    FPHIX = construct_FPHIX(x, y1, y2, y3, hy, hz, kappa, Nx, Ny, Nz);
    FPHIY = construct_FPHIY(x, y1, y2, y3, hy, hz, kappa, Nx, Ny, Nz);
    FPHIZ = construct_FPHIZ(x, y1, y2, y3, hy, hz, kappa, Nx, Ny, Nz);
    
    % create column vectors
    u_x = cube2column(x);
    u_y1 = cube2column(y1);
    u_y2 = cube2column(y2);
    u_y3 = cube2column(y3);
    
    % remove boundary rows (zeros) - NO EQUATIONS AT BOUNDARY
    zero_rows = ~any(LPSIX,2);
    LPSIX(zero_rows, :) = [];
    LPSIY(zero_rows, :) = [];
    LPSIZ(zero_rows, :) = [];

    LPHIX(zero_rows, :) = [];
    LPHIY(zero_rows, :) = [];
    LPHIZ(zero_rows, :) = [];

    FPSI(zero_rows, :) = [];
    FPHIX(zero_rows, :) = [];
    FPHIY(zero_rows, :) = [];
    FPHIZ(zero_rows, :) = [];

    
    %
    dPsidt = D*(LPSIX/hx^2 + LPSIY/hy^2 + LPSIZ/hz^2)*u_x + FPSI;
    dPhidtX = D*(0 + LPHIY/hy^2 + LPHIZ/hz^2)*u_y1 + FPHIX;
    dPhidtY = D*(LPHIX/hx^2 + 0 + LPHIZ/hz^2)*u_y2 + FPHIY;
    dPhidtZ = D*(LPHIX/hx^2 + LPHIY/hy^2 + 0)*u_y3 + FPHIZ;
    


    F = [dPsidt; dPhidtX; dPhidtY; dPhidtZ];
end