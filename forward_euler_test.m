
clear all
close all
eval_u = "analytical_u_xyz";

eval_f = "eval_f";


p.kappa=5;
p.Nx = 100;
p.Ny = 100;
p.Nz = 10;
p.hx = 1e-2;
p.hy = 1e-2;
p.hz = 1e-2;


h = 1;
for k = 2 : p.Nz
    for j = 2 : p.Ny
        for i = 2 : p.Nx
            p.M2(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.m2(h) = i-1 + (p.Nx-1)*(j-2)+(p.Nx-1)*(p.Ny-1)*(k-2);
            h = h + 1;
        end
    end
end
% p.M2 = sub2ind([(p.Nx+1) (p.Nx+1) (p.Ny+1)], [2:p.Nx], [2:p.Ny], [2:p.Nz]);
% p.m2 = sub2ind([(p.Nx-1) (p.Nx-1) (p.Ny-1)], [1:p.Nx-1], [1:p.Ny-1], [1:p.Nz-1]);

h = 1;
for k = 1 : p.Nz+1
    for j = 1 : p.Ny+1
        for i = 1 : p.Nx+1
            p.m(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1); %this is just 1:N^3
            h = h + 1;
        end
    end
end
% 
% p.M2 = sub2ind([(p.Nx+1) (p.Nx+1) (p.Ny+1)], [1:p.Nx+1], [1:p.Ny+1], [1:p.Nz+1]);

% x
h = 1;
h_int = 1;
for k = 1 : p.Nz+1
    for j = 1 : p.Ny+1
        p.m1x(h) = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.m2x(h) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.mNx(h) = p.Nx + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.mNxp1(h) = p.Nx+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        h = h + 1;
        if k ~= 1 && k ~= p.Nz+1 && j ~= 1 && j ~= p.Ny+1
            p.m1x_int(h_int) = 1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.m2x_int(h_int) = 2 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.mNx_int(h_int) = p.Nx + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.mNxp1_int(h_int) = p.Nx+1 + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            h_int = h_int + 1;
        end           
    end
end


% y
h = 1;
h_int = 1;
for k = 1 : p.Nz+1
    for i = 1 : p.Nx+1 
        p.m1y(h) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
        p.m2y(h) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.mNy(h) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
        p.mNyp1(h) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
        h = h + 1;
        if k ~= 1 && k ~= p.Nz+1 && i ~= 1 && i ~= p.Nx+1
            p.m1y_int(h_int) = i + (p.Nx+1)*(p.Ny+1)*(k-1);
            p.m2y_int(h_int) = i + (p.Nx+1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.mNy_int(h_int) = i + (p.Nx+1)*(p.Ny-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            p.mNyp1_int(h_int) = i + (p.Nx+1)*p.Ny+(p.Nx+1)*(p.Ny+1)*(k-1);
            h_int = h_int + 1;
        end           
    end
end

% z
h = 1;
h_int = 1;
for j = 1 : p.Ny+1
    for i = 1 : p.Nx+1
        p.m1z(h) = i + (p.Nx+1)*(j-1);
        p.m2z(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
        p.mNz(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
        p.mNzp1(h) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*p.Nz;
        h = h + 1;
        if j ~= 1 && j ~= p.Ny+1 && i ~= 1 && i ~= p.Nx+1
            p.m1z_int(h_int) = i + (p.Nx+1)*(j-1);
            p.m2z_int(h_int) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1);
            p.mNz_int(h_int) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(p.Nz-1);
            p.mNzp1_int(h_int) = i + (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*p.Nz;
            h_int = h_int + 1;
        end
    end
end

x = sparse(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),1, sqrt(1/((p.Nx-1)*(p.Ny-1)*(p.Nz-1))));
y1 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y2 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);
y3 = sparse((p.Nx-1)*(p.Ny-1)*(p.Nz-1),1);

x_start = [x;y1;y2;y3];

t_start=0;
t_stop=0.001;
max_dt_FE = .0005;

[X] = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,max_dt_FE,1);


%% 
visualizeNetwork(X,p)
