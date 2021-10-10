Nx = 4;
Ny = 4;
Nz = 4;
x = zeros(Nx, Ny, Nz);
y1 = ones(Nx, Ny, Nz);
y2 = ones(Nx, Ny, Nz);
y3 = ones(Nx, Ny, Nz);
hx = 1;
hy = 1;
hz = 1;
kappa=1;

LPHIX = construct_LPHIX(hx,hy,hz,kappa,Nx,Ny,Nz);
LPHIY = construct_LPHIY(hx,hy,hz,kappa,Nx,Ny,Nz);
LPHIZ = construct_LPHIZ(hx,hy,hz,kappa,Nx,Ny,Nz);

LPSIX = construct_LPSIX(y1,Nx, Ny, Nz)./hx^2;
LPSIY = construct_LPSIY(y2,Nx, Ny, Nz)./hy^2;
LPSIZ = construct_LPSIZ(y3,Nx, Ny, Nz)./hz^2;

FPSI = construct_FPSI(x, Nx, Ny, Nz);
D = kappa^2;

dPsi = D*(LPSIX+LPSIY+LPSIZ) + FPSI

figure(1)
imagesc(abs(dPsi))
