function [u,p] = eval_u(t,X,p)

[Bx0, By0, Bz0] = eval_Bfield(X,p);


Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);

classicBx = [p.m1y_int p.m1z_int p.mNyp1_int, p.mNzp1_int];
classicBy = [p.m1z_int p.m1x_int p.mNzp1_int, p.mNxp1_int];
classicBz = [p.m1x_int p.m1y_int p.mNxp1_int, p.mNyp1_int];

if t > 0 && t <= p.t_stop
    p.appliedBx = 0; 
    p.appliedBy = 0; 
    p.appliedBz = p.magBz; 

% elseif t > p.t_stop/2 && t <= p.t_stop
%     p.appliedBx = p.magBx; 
%     p.appliedBy = 0; 
%     p.appliedBz = p.magBz; 
%     
% elseif t > p.t_stop/2 && t <= p.t_stop*3/4
%     p.appliedBx = 0; 
%     p.appliedBy = 0; 
%     p.appliedBz = 1/(p.hx*p.hy)*p.magBz*2; 
    
else
    p.appliedBx = 0;
    p.appliedBy = 0;
    p.appliedBz = 0;

end



%% normalization

%%%%%%% x
if p.appliedBx>0
    sumBx = (2*(p.Nx-2)*(p.Ny-2) + 2*(p.Nx-2)*(p.Nz-2))*p.appliedBx - sum(Bx0);
    sumBx0 = (2*(p.Nx-2)*(p.Ny-2) + 2*(p.Nx-2)*(p.Nz-2))*p.appliedBx;
    p.Brealx = p.appliedBx*abs(sumBx)/abs(sumBx0);
    Bx = Bx + sparse(classicBx,1,p.Brealx,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else
    p.Brealx = 0;
    Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
end


%%%%%%% y
if p.appliedBy>0
    sumBy = (2*(p.Ny-2)*(p.Nx-2) + 2*(p.Ny-2)*(p.Nz-2))*p.appliedBy - sum(By0);
    sumBy0 = (2*(p.Ny-2)*(p.Nx-2) + 2*(p.Ny-2)*(p.Nz-2))*p.appliedBy;
    p.Brealy = p.appliedBy*abs(sumBy)/abs(sumBy0);
    By = By + sparse(classicBy,1,p.Brealy,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else
    p.Brealy = 0;
    By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
end

if p.appliedBz>0
    sumBz = (2*(p.Nz-2)*(p.Nx-2) + 2*(p.Nz-2)*(p.Ny-2))*p.appliedBz - sum(Bz0);
    sumBz0 = (2*(p.Nz-2)*(p.Nx-2) + 2*(p.Nz-2)*(p.Ny-2))*p.appliedBz;
    p.Brealz = p.appliedBz*abs(sumBz)/abs(sumBz0);
    Bz = Bz + sparse(classicBz,1,p.Brealz,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else
    p.Brealz = 0;
    Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
end

u.Bx = Bx;
u.By = By;
u.Bz = Bz;

p.u = u;
end
