function p = contruct_indices(p)

    mk = (p.Nx+1)*(p.Ny+1);
    mj = (p.Nx+1);

    p.m2 = 1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1);
    p.m = 1:(p.Nx+1)*(p.Ny+1)*(p.Nz+1);
   
    
    % indices for vectors
    h_x = 1;
    h_int_x = 1;
    h_y = 1;
    h_int_y = 1;
    h_z = 1;
    h_int_z = 1;
    h_M2 = 1;

    for k = 1 : p.Nz+1
        for j = 1 : p.Ny+1
            for i = 1 : p.Nx+1

                if k ~= 1 && k ~= p.Nz+1 && j ~= 1 && j ~= p.Ny+1 && i ~= 1 && i ~= p.Nx+1
                    p.M2(h_M2) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    h_M2 = h_M2 + 1;
                end

                if i == 1
                    % x
                    p.m1x(h_x) = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    p.m2x(h_x) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    p.mNx(h_x) = p.Nx + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    p.mNxp1(h_x) = p.Nx+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    h_x = h_x + 1;
                    if k ~= 1 && k ~= p.Nz+1 && j ~= 1 && j ~= p.Ny+1
                        p.m1x_int(h_int_x) = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.m2x_int(h_int_x) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNx_int(h_int_x) = p.Nx + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNxp1_int(h_int_x) = p.Nx+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        h_int_x = h_int_x + 1;
                    end 
                end

                if j == 1
                    % y
                    p.m1y(h_y) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
                    p.m2y(h_y) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    p.mNy(h_y) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    p.mNyp1(h_y) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
                    h_y = h_y + 1;
                    if k ~= 1 && k ~= p.Nz+1 && i ~= 1 && i ~= p.Nx+1
                        p.m1y_int(h_int_y) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
                        p.m2y_int(h_int_y) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNy_int(h_int_y) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNyp1_int(h_int_y) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
                        h_int_y = h_int_y + 1;
                    end           
                end

                if k == 1
                    % z
                    p.m1z(h_z) = i + (p.Nx+1)*(j-1);
                    p.m2z(h_z) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
                    p.mNz(h_z) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
                    p.mNzp1(h_z) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*p.Nz;
                    h_z = h_z + 1;
                    if j ~= 1 && j ~= p.Ny+1 && i ~= 1 && i ~= p.Nx+1
                        p.m1z_int(h_int_z) = i + (p.Nx+1)*(j-1);
                        p.m2z_int(h_int_z) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
                        p.mNz_int(h_int_z) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
                        p.mNzp1_int(h_int_z) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*p.Nz;
                        h_int_z = h_int_z + 1;
                    end
                end

            end
        end
    end
    
    
%     % indices for L
%     m = p.M2; 
%     dim_L = [(p.Nx+1)*(p.Ny+1)*(p.Nz+1),(p.Nx+1)*(p.Ny+1)*(p.Nz+1)];
%     p.L_m = sub2ind(dim_L,m,m);
%     p.L_pmj = sub2ind(dim_L,m,m+mj);
%     p.L_mmj = sub2ind(dim_L,m,m-mj);
%     p.L_pmk = sub2ind(dim_L,m,m+mk);
%     p.L_mmk = sub2ind(dim_L,m,m-mk);
%     p.L_p1 = sub2ind(dim_L,m,m+1);
%     p.L_m1 = sub2ind(dim_L,m,m-1);
end

