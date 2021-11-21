function [u,p] = eval_u(t,X,p)

[Bx0, By0, Bz0] = eval_Bfield(X,p);

Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);



classicBx = [p.m1y_int p.m1z_int p.mNyp1_int, p.mNzp1_int];
classicBy = [p.m1z_int p.m1x_int p.mNzp1_int, p.mNxp1_int];
classicBz = [p.m1x_int p.m1y_int p.mNxp1_int, p.mNyp1_int];

if t > 0 && t < 1
    Bx(p.M2B) = Bx0;
    By(p.M2B) = By0;
    Bz(p.M2B) = Bz0;
    p.appliedBz = 1/(p.hx*p.hy)*p.magBz*t;
    Bx = Bx + sparse(classicBx,1,1/(p.hy*p.hz)*p.magBx*t,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    By = By + sparse(classicBy,1,1/(p.hz*p.hx)*p.magBy*t,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    Bz = Bz + sparse(classicBz,1,p.appliedBz,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
elseif t > 1
    Bx(p.M2B) = Bx0;
    By(p.M2B) = By0;
    Bz(p.M2B) = Bz0;
    p.appliedBz = 0.03;
    Bx = Bx + sparse(classicBx,1,1/(p.hy*p.hz)*p.magBx/10,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    By = By + sparse(classicBy,1,1/(p.hz*p.hx)*p.magBy/10,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    Bz = Bz + sparse(classicBz,1,p.appliedBz,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else

    Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    Bx(p.M2B) = Bx(p.M2B) + Bx0;
    By(p.M2B) = By(p.M2B) + By0;
    Bz(p.M2B) = Bx(p.M2B) + Bz0;
    p.appliedBz = 0;
end
u.Bx = Bx;
u.By = By;
u.Bz = Bz;

p.u = u;
end

% 
% function u = eval_u(t,X,p)
% 
% [Bx0, By0, Bz0] = eval_Bfield(X,p);
% 
% Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
% By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
% Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
% 
% 
% 
% classicBx = [p.m1y_int p.m1z_int p.mNyp1_int, p.mNzp1_int];
% classicBy = [p.m1z_int p.m1x_int p.mNzp1_int, p.mNxp1_int];
% classicBz = [p.m1x_int p.m1y_int p.mNxp1_int, p.mNyp1_int];
% 
% if t > 0 && t < 1
%     Bx(p.M2B) = Bx0;
%     By(p.M2B) = By0;
%     Bz(p.M2B) = Bz0;
% 
%     Bx = Bx + sparse(classicBx,1,1/(p.hy*p.hz)*p.magBx*t,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     By = By + sparse(classicBy,1,1/(p.hz*p.hx)*p.magBy*t,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     Bz = Bz + sparse(classicBz,1,1/(p.hx*p.hy)*p.magBz*t,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
% else
%     Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
% end
% max(Bz)
% u.Bx = Bx;
% u.By = By;
% u.Bz = Bz;
% 
% p.u = u;
% end