function u = eval_u(t,p)

x = sparse(M2, 1, x_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
y1 = sparse(M2, 1, y1_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
y2 = sparse(M2, 1, y2_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);
y3 = sparse(M2, 1, y3_int, (p.Nx+1)*(p.Ny+1)*(p.Nz+1), 1);

if t <0
   u.Bx = 0;
   u.By = 0;
   u.Bz = 0;
else 
   u.Bx = 0;
   u.By = 0;
   u.Bz = 1;
end

end