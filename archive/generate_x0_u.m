function [x0, y1_b, y2_b, y3_b, y1_0, y2_0, y3_0] = generate_x0_u(Bx, Nx, Ny, Nz, hy, hz)
    row_B = [1 : Bx*hy*hz : Bx*hy*hz*Ny];
    y3_b = row_B'*ones(1,Nz-1);
    x0 = ones(Nx, Ny, Nz);
    y1_0 = ones(Nx-1, Ny, Nz);
    y2_0 = ones(Nx, Ny-1, Nz);
    y3_0 = ones(Nx, Ny, Nz-1);

    % phi_z in first layer yz
    y3_0(1, :, :) = y3_b;
    y3_0(2, :, :) = y3_b;

    % phi_z in last layer yz
    y3_0(Nx-1, :, :) = y3_b;
    y3_0(Nx, :, :) = y3_b;
end

