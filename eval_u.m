function u = eval_u(t,p)

Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);


classicBx = [p.m1y_int p.m1z_int p.mNyp1_int, p.mNzp1_int];
classicBy = [p.m1z_int p.m1x_int p.mNzp1_int, p.mNxp1_int];
classicBz = [p.m1x_int p.m1y_int p.mNxp1_int, p.mNyp1_int];
if t <0
    Bx = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    By = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    Bz = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
else 
    Bx = Bx + sparse(classicBx,1,p.magBx,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    By = By + sparse(classicBy,1,p.magBy,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
    Bz = Bz + sparse(classicBz,1,p.magBz,(p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
end

u.Bx = Bx;
u.By = By;
u.Bz = Bz;

p.u = u;
end