
% Project Milestone 1

hx = .2;
hy = .2;
% hz = .1;
Nx = 1/hx;
Ny = 1/hy;
K = 1;
Ux = sym('Ux', [Nx;1]);
Uy = sym('Uy', [Ny;1]);

PSIx = sym('PSIx', [Nx;1]);
PSIy = sym('PSIy', [Nx;1]);

LX_phi = eval_L(Ux, hx);
LY_phi = eval_L(Uy, hx);



% abs(LX_phi)

dPSI = LX_phi*PSIx + LY_phi*PSIy