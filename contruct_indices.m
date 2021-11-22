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
    h_M2B = 1;
    
    h_edge_x = 1;
    h_edge_y = 1;
    h_edge_z = 1;
    
    h_edge_x21 = 1;
    h_edge_y21 = 1;
    h_edge_z21 = 1;
    
    h_edge_x22 = 1;
    h_edge_y22 = 1;
    h_edge_z22 = 1;

    h_checkx1 = 1;
    h_checkx2 = 1;
    h_checky1 = 1;
    h_checky2 = 1;
    h_checkz1 = 1;
    h_checkz2 = 1;
    
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
                    p.mNx(h_x) = p.Ny + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    p.mNxp1(h_x) = p.Ny+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    h_x = h_x + 1;
                    if k ~= 1 && k ~= p.Nz+1 && j ~= 1 && j ~= p.Nx+1
                        p.m1x_int(h_int_x) = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.m2x_int(h_int_x) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNx_int(h_int_x) = p.Ny + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNxp1_int(h_int_x) = p.Ny+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        if mod(j,2)==0
                            p.m2x_int_cb1(h_checkx1) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                            p.mNx_int_cb1(h_checkx1) = p.Ny + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                            h_checkx1 = h_checkx1+1;
                        else
                            p.m2x_int_cb2(h_checkx2) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                            p.mNx_int_cb2(h_checkx2) = p.Ny + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);                
                            h_checkx2 = h_checkx2+1;
                        end
                        h_int_x = h_int_x + 1;
                    end 
                    
                    if (j==1 || j==p.Ny+1) 
                        p.m1xedge(h_edge_x) = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNxedgep1(h_edge_x) =  p.Ny+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        h_edge_x = h_edge_x+1;
                    end
                    

                    if (j==2) && (k ~= 1 && k ~= p.Nz+1)
                        p.m2xedge1(h_edge_x21) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNxedge1(h_edge_x21) =  p.Ny + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        h_edge_x21 = h_edge_x21+1;
                    end     
                    
                    if (j==p.Ny) && (k ~= 1 && k ~= p.Nz+1)
                        p.m2xedge2(h_edge_x22) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNxedge2(h_edge_x22) =  p.Ny + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        h_edge_x22 = h_edge_x22+1;
                    end
                end

                if j == 1
                    % y
                    p.m1y(h_y) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
                    p.m2y(h_y) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    p.mNy(h_y) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                    p.mNyp1(h_y) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
                    h_y = h_y + 1;
                    if k ~= 1 && k ~= p.Nz+1 && i ~= 1 && i ~= p.Ny+1
                        p.m1y_int(h_int_y) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
                        p.m2y_int(h_int_y) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNy_int(h_int_y) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNyp1_int(h_int_y) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
                        
                        if mod(k,2)==0
                            p.m2y_int_cb1(h_checky1) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                            p.mNy_int_cb1(h_checky1) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                            h_checky1 = h_checky1+1;
                        else
                            p.m2y_int_cb2(h_checky2) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                            p.mNy_int_cb2(h_checky2) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);                    
                            h_checky2 = h_checky2+1;
                        end
                        
                        h_int_y = h_int_y + 1;
                    end     
                    
                    if (k==1 || k==p.Nz+1)
                        p.m1yedge(h_edge_y) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNyedgep1(h_edge_y) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
                        h_edge_y = h_edge_y+1;
                    end
                    
                    if (k==2) && (i ~= 1 && i ~= p.Ny+1)
                        p.m2yedge1(h_edge_y21) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNyedge1(h_edge_y21) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        h_edge_y21 = h_edge_y21+1;
                    end
                    if (k==p.Nz) && (i ~= 1 && i ~= p.Ny+1)
                        p.m2yedge2(h_edge_y22) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        p.mNyedge2(h_edge_y22) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
                        h_edge_y22 = h_edge_y22+1;
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
                        
                        if mod(i,2)==0
                            p.m2z_int_cb1(h_checkz1) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
                            p.mNz_int_cb1(h_checkz1) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
                            h_checkz1 = h_checkz1+1;
                        else
                            p.m2z_int_cb2(h_checkz2) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
                            p.mNz_int_cb2(h_checkz2) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);                        
                            h_checkz2 = h_checkz2+1;
                        end
                        h_int_z = h_int_z + 1;
                    end
                                        
                    if (i==1 || i==p.Nz+1)
                        p.m1zedge(h_edge_z) = i + (p.Nx+1)*(j-1);
                        p.mNzedgep1(h_edge_z) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*p.Nz;
                        h_edge_z = h_edge_z+1;
                    end
                    if (i==2) && (j ~= 1 && j ~= p.Ny+1)
                        p.m2zedge1(h_edge_z21) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
                        p.mNzedge1(h_edge_z21) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
                        h_edge_z21 = h_edge_z21+1;
                    end
                    if (i==p.Nz) && (j ~= 1 && j ~= p.Ny+1)
                        p.m2zedge2(h_edge_z22) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
                        p.mNzedge2(h_edge_z22) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
                        h_edge_z22 = h_edge_z22+1;
                    end
                end

            end
        end
    end
    
    for k = 1 : p.Nz
        for j = 1 : p.Ny
            for i = 1 : p.Nx

%                 if k ~= 1 && k ~= p.Nz && j ~= 1 && j ~= p.Ny && i ~= 1 && i ~= p.Nx
                if k < p.Nz-1 && j < p.Ny-1 && i < p.Nx-1
                    p.M2B(h_M2B) = i + (p.Nx-1)*(j-1)+(p.Nx-1)*(p.Ny-1)*(k-1);
                    h_M2B = h_M2B + 1;
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

