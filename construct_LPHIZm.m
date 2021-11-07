function LPHI = construct_LPHIZm(p)
hx = p.hx;
hy = p.hy;
% hz = p.hz; 
kappa = p.kappa; 
% Nx = p.Nx;
% Ny = p.Ny;
% Nz = p.Nz;
LPHI = sparse((p.Nx+1)*(p.Ny+1)*(p.Nz+1), (p.Nx+1)*(p.Ny+1)*(p.Nz+1));
% 
% mk = (p.Nx+1)*(p.Ny+1);
% mj = (p.Nx+1);
% m = p.M2;
%    
%             
            
% Diagonal
LPHI(p.L_m) = - 2*(kappa^2/hx^2 + kappa^2/hy^2);

% i TERMS
if Nx > 1
    LPHI(p.L_p1) = kappa^2/hx^2;
    LPHI(p.L_m1) = kappa^2/hx^2;
end

            
% k TERMS
if Ny > 1
    LPHI(p.L_pmj) = kappa^2/hy^2;
    LPHI(p.L_mmj) = kappa^2/hy^2;
end


end