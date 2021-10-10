% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)
% 
% x(i) = x(i+1)*exp(-1i*x(N+i)); %psi1
% x(N) = x(i-1)*exp(1i*x(N+i-1)); %psiN

function [FPHI_x, FPHI_y, FPHI_z] = construct_FPHI(x, y1, y2, y3, hx, hy, hz, kappa, Nx, Ny, Nz)

    for k = 2 : Nz-1
        for j = 2 : Ny-1
            for i = 2 : Nx-1
                                           
                FPHI_x(i,j,k) = (kappa^2/hy^2)*(-y2(i+1,j,k)+y2(i,j,k)+y2(i+1,j-1,k)-y2(i,j-1,k))...
                                +(kappa^2/hz^2)*(-y3(i+1,j,k)+y3(i,j,k)+y3(i+1,j,k-1)-y3(i,j,k-1))...
                                +imag(exp(-1i*y1(i,j,k))*conj(x(i,j,k))*x(i+1,j,k));

                FPHI_y(i,j,k) = (kappa^2/hz^2)*(-y3(i,j+1,k)+y3(i,j,k)+y3(i,j+1,k-1)-y3(i,j,k-1))...
                                +(kappa^2/hx^2)*(-y1(i,j+1,k)+y1(i,j,k)+y1(i-1,j+1,k)-y1(i-1,j,k))...
                                +imag(exp(-1i*y2(i,j,k))*conj(x(i,j,k))*x(i,j+1,k));

                FPHI_z(i,j,k) =(kappa^2/hx^2)*(-y1(i,j,k+1)+y1(i,j,k)+y1(i-1,j,k+1)-y1(i-1,j,k))...
                                +(kappa^2/hy^2)*(-y2(i,j,k+1)+y2(i,j,k)+y2(i,j-1,k+1)-y2(i,j-1,k))...
                                +imag(exp(-1i*y3(i,j,k))*conj(x(i,j,k))*x(i,j,k+1));          
            
            end
        end
    end
end                            
