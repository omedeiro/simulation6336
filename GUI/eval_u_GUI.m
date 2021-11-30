function [u,p] = eval_u_GUI(t,X,p,n)

% global cord
% global click
% global click_location
% x = round(cord(1,1),0);
% y = round(cord(1,2),0);
% int_val = y + (p.Nx-1)*(x-1);
% 
% if click == 1
%     click_location = [click_location int_val];
%     disp("x: "+x + " y: " + y + " loc: " + click_location(end) + " click = " + click)
% end
[Bx0, By0, Bz0] = eval_Bfield(X,p);

p.BXT(:,n) = Bx0;
p.BYT(:,n) = By0;
p.BZT(:,n) = Bz0;

Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);

classicBx = [p.m1y_int p.m1z_int p.mNyp1_int, p.mNzp1_int];
classicBy = [p.m1z_int p.m1x_int p.mNzp1_int, p.mNxp1_int];
classicBz = [p.m1x_int p.m1y_int p.mNxp1_int, p.mNyp1_int];

if t > 0 && t <= p.t_stop
    p.appliedBx = p.magBx; 
    p.appliedBy = p.magBy; 
    p.appliedBz = p.magBz; 

% elseif t > p.t_stop/2 && t <= p.t_stop*3/4
%     p.appliedBx = 0; 
%     p.appliedBy = 0; 
%     p.appliedBz = 0;     
    
else
    p.appliedBx = 0;
    p.appliedBy = 0;
    p.appliedBz = 0;

end



%% normalization

%%%%%%% x
if abs(p.appliedBx)>0
    sumBx = (2*(p.Nx-2)*(p.Ny-2) + 2*(p.Nx-2)*(p.Nz-2))*p.appliedBx - sum(Bx0);
    sumBx0 = (2*(p.Nx-2)*(p.Ny-2) + 2*(p.Nx-2)*(p.Nz-2))*p.appliedBx;
    p.Brealx = p.appliedBx;%*abs(sumBx)/abs(sumBx0);
    Bx = Bx + sparse(classicBx,1,p.Brealx,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     Bx = Bx + sparse(p.m1z_int(click_location),1,p.Brealx,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     Bx = Bx + sparse(p.mNzp1_int(click_location),1,p.Brealx,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else
    p.Brealx = 0;
    Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
end


%%%%%%% y
if abs(p.appliedBy)>0
    sumBy = (2*(p.Ny-2)*(p.Nx-2) + 2*(p.Ny-2)*(p.Nz-2))*p.appliedBy - sum(By0);
    sumBy0 = (2*(p.Ny-2)*(p.Nx-2) + 2*(p.Ny-2)*(p.Nz-2))*p.appliedBy;
    p.Brealy = p.appliedBy;%*abs(sumBy)/abs(sumBy0);
    By = By + sparse(classicBy,1,p.Brealy,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     By = By + sparse(p.m1z_int(click_location),1,p.Brealy,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     By = By + sparse(p.mNzp1_int(click_location),1,p.Brealy,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else
    p.Brealy = 0;
    By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
end

if abs(p.appliedBz)>0
    sumBz = (2*(p.Nz-2)*(p.Ny-2) + 2*(p.Nz-2)*(p.Nx-2))*p.appliedBz - sum(Bz0);
    sumBz0 = (2*(p.Nz-2)*(p.Ny-2) + 2*(p.Nz-2)*(p.Nx-2))*p.appliedBz;
    p.Brealz = p.appliedBz;%*abs(sumBz)/abs(sumBz0);
%     p.Brealz = p.appliedBz*abs(sumBz)/abs(sumBz0);
%     p.Brealz = real(p.appliedBz - max(Bz0));
    Bz = Bz + sparse(classicBz,1,p.Brealz,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     Bz = Bz + sparse(p.m1z_int(click_location),1,p.Brealz,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
%     Bz = Bz + sparse(p.mNzp1_int(click_location),1,p.Brealz,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else
    p.Brealz = 0;
    Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
end

u.Bx = Bx;
u.By = By;
u.Bz = Bz;

p.u = u;
end
