% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)
% 
% x(i) = x(i+1)*exp(-1i*x(N+i)); %psi1
% x(N) = x(i-1)*exp(1i*x(N+i-1)); %psiN

% function F = analytical_f_xyz(x_int, y1_int, y2_int, y3_int, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
% function F = analytical_f_xyz(X, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
function F = eval_f(X, p, u)

    colN = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
    x = X(1:colN);
    y1 = X(colN+1:2*colN);
    y2 = X(2*colN+1:3*colN);
    y3 = X(3*colN+1:4*colN);


    %%%%%%%% _int has size (p.Nx-1, p.Ny-1, p.Nz-1)
    x = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    y1 = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    y2 = zeros(p.Nx+1, p.Ny+1, p.Nz+1); 
    y3 = zeros(p.Nx+1, p.Ny+1, p.Nz+1);
    
    x(2:p.Nx, 2:p.Ny, 2:p.Nz) = x_int;
    y1(2:p.Nx, 2:p.Ny, 2:p.Nz) = y1_int;
    y2(2:p.Nx, 2:p.Ny, 2:p.Nz) = y2_int;
    y3(2:p.Nx, 2:p.Ny, 2:p.Nz) = y3_int;
    
    
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
        y1(p.Nx+1, :, :) = ones(p.Ny+1, p.Nz+1);

        %%%% phi_z in last layer yz
        % first face
        y3(p.Nx+1, :, :) = y3_b;
        y2(p.Nx+1, :, :) = ones(p.Ny+1, p.Nz+1);
    %     % second face
    %     y3(p.Nx+1, :, :) = y3_b;
    %     y2(p.Nx+1, :, :) = zeros(p.Ny+1, p.Nz+1);


    else
        % Magnetic field boundary conditions eq 37

        for kk = 2:p.Nz
            for ii = 2:p.Nx
%                 m = Nx*(j-1)+Nx*Ny*(k-1);

                y1(ii,1,kk) = u.Bz*p.hx*p.hy + y1(ii,2,kk) + y2(ii,1,kk) - y2(ii+1,1,kk);
                y1(ii,p.Ny+1,kk) = u.Bz*p.hx*p.hy + y1(ii,p.Ny,kk) - y2(ii,p.Ny+1,kk) + y2(ii+1,p.Ny+1,kk);
                
                y3(ii,1,kk) = -u.Bx*p.hy*p.hz + y2(ii,1,kk) - y2(ii,1,kk+1) + y3(ii,2,kk);
                y3(ii,p.Ny+1,kk) = -u.Bx*p.hy*p.hz + y2(ii,p.Ny+1,kk) - y2(ii,p.Ny+1,kk+1) + y3(ii,p.Ny,kk);
                
            end
        end
        for jj = 2:p.Ny
            for ii = 2:p.Nx
                
                y2(ii,jj,1) = u.Bx*p.hy*p.hz + y2(ii,jj,2) + y3(ii,jj,1) - y3(ii,jj+1,1);
                y2(ii,jj,p.Nz+1) = u.Bx*p.hy*p.hz + y2(ii,jj,p.Nz) + y3(ii,jj,p.Nz+1) - y3(ii,jj+1,p.Nz+1);
                
                y1(ii,jj,1) = -u.By*p.hz*p.hx + y3(ii,jj,1) - y3(ii+1,jj,1) + y1(ii,jj,2);
                y1(ii,jj,p.Nz+1) = -u.By*p.hz*p.hx + y3(ii,jj,p.Nz+1) - y3(ii+1,jj,p.Nz+1) + y1(ii,jj,p.Nz);
          
            end
        end
        for kk = 2:p.Nz
            for jj = 2:p.Ny              
                
                y2(1,jj,kk) = -u.Bz*p.hx*p.hy + y1(1,jj,kk) - y1(1,jj+1,kk) + y2(2,jj,kk);
                y2(p.Nx+1,jj,kk) = -u.Bz*p.hx*p.hy + y1(p.Nx+1,jj,kk) - y1(p.Nx+1,jj+1,kk) + y2(p.Nx,jj,kk);
                
                y3(1,jj,kk) = u.By*p.hz*p.hx + y3(2,jj,kk) +y1(1,jj,kk) - y1(1,jj,kk+1);
                y3(p.Nx+1,jj,kk) = u.By*p.hz*p.hx + y3(p.Nx,jj,kk) +y1(p.Nx+1,jj,kk) - y1(p.Nx+1,jj,kk+1);
                
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
    
    % create column vectors
    u_x = sparse(cube2column(x));
    u_y1 = sparse(cube2column(y1));
    u_y2 = sparse(cube2column(y2));
    u_y3 = sparse(cube2column(y3));
    
    % remove boundary rows (zeros) - NO EQUATIONS AT BOUNDARY
    
    LPSIX(~any(LPSIX,2), :) = [];
    LPSIY(~any(LPSIY,2), :) = [];
    LPSIZ(~any(LPSIZ,2), :) = [];

    LPHIX(~any(LPHIX,2), :) = [];
    LPHIY(~any(LPHIY,2), :) = [];
    LPHIZ(~any(LPHIZ,2), :) = [];

    %
    dPsidt = D*(LPSIX/p.hx^2 + LPSIY/p.hy^2 + LPSIZ/p.hz^2)*u_x + FPSI;
    dPhidtX = D*(LPHIY./p.hy^2 + LPHIZ./p.hz^2)*u_y1 + FPHIX;
    dPhidtY = D*(LPHIX./p.hx^2 + LPHIZ./p.hz^2)*u_y2 + FPHIY;
    dPhidtZ = D*(LPHIX./p.hx^2 + LPHIY./p.hy^2)*u_y3 + FPHIZ;
    


    F = [dPsidt; dPhidtX; dPhidtY; dPhidtZ];
%     F = abs(F);
end