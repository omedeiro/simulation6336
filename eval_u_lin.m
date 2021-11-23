function [u,p] = eval_u(t,X,p)

[Bx0, By0, Bz0] = eval_Bfield(X,p);


Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);

% classicBx = [p.m1y_int p.m1z_int p.mNyp1_int, p.mNzp1_int];
% classicBy = [p.m1z_int p.m1x_int p.mNzp1_int, p.mNxp1_int];
classicBz = [p.m1x_int p.m1y_int p.mNxp1_int, p.mNyp1_int];

if t > 0 && t < 50
    p.appliedBz = 1/(p.hx*p.hy)*p.magBz*t; 
elseif t > 50 && t < 100
    p.appliedBz =  1/(p.hx*p.hy)*(p.magBz*(100-t)); 
% elseif t > 300 && t < 310
%     p.appliedBz = 1/(p.hx*p.hy)*p.magBz*(t-290); 
% elseif t > 310 && t < 600
%     p.appliedBz =  1/(p.hx*p.hy)*p.magBz*20; 
% elseif t > 600 && t < 610
%     p.appliedBz =  1/(p.hx*p.hy)*p.magBz*(t-580); 
% elseif t > 610
%     p.appliedBz =  1/(p.hx*p.hy)*p.magBz*30; 
   
else
    p.appliedBz = 0;
end

if p.appliedBz>0
    sumBz = (2*(p.Nz+1)*(p.Nx+1) + 2*(p.Nz+1)*(p.Ny+1))*p.appliedBz - sum(Bz0);
    sumBz0 = (2*(p.Nz+1)*(p.Nx+1) + 2*(p.Nz+1)*(p.Ny+1))*p.appliedBz;
    p.Breal = p.appliedBz*abs(sumBz)/abs(sumBz0);
    Bz = Bz + sparse(classicBz,1,p.Breal,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else
    p.Breal = 0;
    Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
end

u.Bx = Bx;
u.By = By;
u.Bz = Bz;

p.u = u;
end
