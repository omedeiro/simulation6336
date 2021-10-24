function u = analytical_u_xyz(t)
u.By = 0;
u.Bz = 0;
if t <0
   u.Bx = 0;
else 
   u.Bx = 1;
end

end