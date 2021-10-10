% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)
% 
% x(i) = x(i+1)*exp(-1i*x(N+i)); %psi1
% x(N) = x(i-1)*exp(1i*x(N+i-1)); %psiN

function F = analytical_f_xyz(x, y1, y2, y3, Bx, hx, hy, hz, kappa, Nx, Ny, Nz)
    dPsidt = zeros(Nx, Ny, Nz);
    dPhidtX = zeros(Nx, Ny, Nz);
    dPhidtY = zeros(Nx, Ny, Nz);
    dPhidtZ = zeros(Nx, Ny, Nz);
    
    % BCs on Bx
    row_B = [1 : Bx*hy*hz : Bx*hy*hz*Ny];
    y3_b = row_B'*ones(1,Nz-1);

    % phi_z in first layer yz
    y3(1, :, :) = y3_b;
    y3(2, :, :) = y3_b;

    % phi_z in last layer yz
    y3(Nx-1, :, :) = y3_b;
    y3(Nx, :, :) = y3_b;

%     %%% BOUNDARY CONDITIONS %%%
%     for k = 1:N
%         for j = 1:N
%             for i = 1:N
%                 x(i,j,1) = x(i,j,N-1); %34
%                 y1(i,j,1) = y1(i,j,N-1); %34
%                 x(i,j,2) = x(i,j,N); %34
%                 y1(i,j,2) = y1(i,j,N); %34
%                 
%                 x(1,j,k) = x(2,j,k)*exp(-1i*y1(1,j,k)); %35
%                 x(N,j,k) = x(N-1,j,k)*exp(1i*y1(N-1,j,k)); %35
%                 x(i,1,k) = x(2,j,k)*exp(-1i*y2(i,1,k));
%                 x(i,N,k) = x(i,N-1,k)*exp(1i*y2(i,N-1,k));
%             end
%         end
%     end
%     
    for k = 1 : Nz-1
        for j = 1 : Ny-1
            for i = 2 : Nx-1
                %%%%%%%% CENTER %%%%%%%%%%%
                dPhidtX(i,j,k) = (kappa^2/hy^2)*(y1(i,j+1,k)-2*y1(i,j,k)+y1(i,j-1,k))...
                                +(kappa^2/hz^2)*(y1(i,j,k+1)-2*y1(i,j,k)+y1(i,j,k-1))...
                                +(kappa^2/hy^2)*(-y2(i+1,j,k)+y2(i,j,k)+y2(i+1,j-1,k)-y2(i,j-1,k))...
                                +(kappa^2/hz^2)*(-y3(i+1,j,k)+y3(i,j,k)+y3(i+1,j,k-1)-y3(i,j,k-1))...
                                +imag(exp(-1i*y1(i,j,k))*conj(x(i,j,k))*x(i+1,j,k));

                dPhidtY(i,j,k) = (kappa^2/hz^2)*(y2(i,j,k+1)-2*y2(i,j,k)+y2(i,j,k-1))...
                                +(kappa^2/hx^2)*(y2(i+1,j,k)-2*y2(i,j,k)+y2(i-1,j,k))...
                                +(kappa^2/hz^2)*(-y3(i,j+1,k)+y3(i,j,k)+y3(i,j+1,k-1)-y3(i,j,k-1))...
                                +(kappa^2/hx^2)*(-y1(i,j+1,k)+y1(i,j,k)+y1(i-1,j+1,k)-y1(i-1,j,k))...
                                +imag(exp(-1i*y2(i,j,k))*conj(x(i,j,k))*x(i,j+1,k));

                dPhidtZ(i,j,k) = (kappa^2/hx^2)*(y3(i+1,j,k)-2*y3(i,j,k)+y3(i-1,j,k))...
                                +(kappa^2/hy^2)*(y3(i,j+1,k)-2*y3(i,j,k)+y3(i,j-1,k))...
                                +(kappa^2/hx^2)*(-y1(i,j,k+1)+y1(i,j,k)+y1(i-1,j,k+1)-y1(i-1,j,k))...
                                +(kappa^2/hy^2)*(-y2(i,j,k+1)+y2(i,j,k)+y2(i,j-1,k+1)-y2(i,j-1,k))...
                                +imag(exp(-1i*y3(i,j,k))*conj(x(i,j,k))*x(i,j,k+1));
            end
        end
    end
    
    F = [dPsidt;dPhidtX;dPhidtY;dPhidtZ];
end                            
%                 %%%%%%%%%%%% ORIGIN %%%%%%%%%%
%                 if i==1 && j==1 && k==1 
%                     dPsidt(i,j,k) = (-x(i,j,k) + exp(-1i*y1(i,j,k))*x(i+1,j,k))/hx^2 ...
%                                     + (-x(i,j,k) + exp(-1i*y2(i,j,k))*x(i,j+1,k))/hy^2 ...
%                                     + (-x(i,j,k) + exp(-1i*y3(i,j,k))*x(i+1,j,k+1))/hz^2 ...
%                                     + (1-x(i,j,k)*conj(x(i,j,k)))*x(i,j,k);
% 
%                     dPhidtX(i,j,k) = (kappa^2/hy^2)*(y1(i,j+1,k)-y1(i,j,k))...
%                                     +(kappa^2/hz^2)*(y1(i,j,k+1)-y1(i,j,k))...
%                                     +(kappa^2/hy^2)*(-y2(i+1,j,k)+y2(i,j,k))...
%                                     +(kappa^2/hz^2)*(-y3(i+1,j,k)+y3(i,j,k))...
%                                     +imag(exp(-1i*y1(i,j,k))*conj(x(i,j,k))*x(i+1,j,k));
% 
%                     dPhidtY(i,j,k) = (kappa^2/hz^2)*(y2(i,j,k+1)-y2(i,j,k))...
%                                     +(kappa^2/hx^2)*(y2(i+1,j,k)-y2(i,j,k))...
%                                     +(kappa^2/hz^2)*(-y3(i,j+1,k)+y3(i,j,k))...
%                                     +(kappa^2/hx^2)*(-y1(i,j+1,k)+y1(i,j,k))...
%                                     +imag(exp(-1i*y2(i,j,k))*conj(x(i,j,k))*x(i,j+1,k));
% 
%                     dPhidtZ(i,j,k) = (kappa^2/hx^2)*(y3(i+1,j,k)-y3(i,j,k))...
%                                     +(kappa^2/hy^2)*(y3(i,j+1,k)-y3(i,j,k))...
%                                     +(kappa^2/hx^2)*(-y1(i,j,k+1)+y1(i,j,k))...
%                                    +(kappa^2/hy^2)*(-y2(i,j,k+1)+y2(i,j,k))...
%                                     +imag(exp(-1i*y3(i,j,k))*conj(x(i,j,k))*x(i,j,k+1));
%                 end
%                 
%                 
%                 
%                 
%                 if i==1 && j==1 && k > 1 && k < N
%                     dPsidt(i,j,k) = (exp(-x(i,j,k) + exp(-1i*y1(i,j,k))*x(i+1,j,k))/hx^2 ...
%                                     +(exp(1i*y2(i,j-1,k))*x(i,j-1,k) - 2*x(i,j,k) + exp(-1i*y2(i,j,k))*x(i,j+1,k))/hy^2 ...
%                                     +(exp(1i*y3(i,j,k-1))*x(i,j,k-1) - 2*x(i,j,k) + exp(-1i*y3(i,j,k))*x(i,j,k+1))/hz^2 ...
%                                     + (1-x(i,j,k)*conj(x(i,j,k)))*x(i,j,k);  
%                                 
%                     dPhidtX(i,j,k) = (kappa^2/hy^2)*(y1(i,j+1,k)-y1(i,j,k))...
%                                     +(kappa^2/hz^2)*(y1(i,j,k+1)-y1(i,j,k))...
%                                     +(kappa^2/hy^2)*(-y2(i+1,j,k)+y2(i,j,k))...
%                                     +(kappa^2/hz^2)*(-y3(i+1,j,k)+y3(i,j,k))...
%                                     +imag(exp(-1i*y1(i,j,k))*conj(x(i,j,k))*x(i+1,j,k));
%                     
%                 end
%                 
%                 if j==1 && k==1
%                     
%                 end
%                 
%                 if k==1 && i==1
%                     
%                 end
%                 
%                 if i==N && j==N && k==N
%                     dPsidt(i,j,k) = (exp(1i*y1(i-1,j,k))*x(i-1,j,k) - x(i,j,k))/hx^2 ...
%                                     +(exp(1i*y2(i,j-1,k))*x(i,j-1,k) - x(i,j,k))/hy^2 ...
%                                     +(exp(1i*y3(i,j,k-1))*x(i,j,k-1) - x(i,j,k))/hz^2 ...
%                                     + (1-x(i,j,k)*conj(x(i,j,k)))*x(i,j,k);  
%                                 
%                     dPhidtX(i,j,k) = (kappa^2/hy^2)*(y1(i,j,k)+y1(i,j-1,k))...
%                                     +(kappa^2/hz^2)*(y1(i,j,k)+y1(i,j,k-1))...
%                                     +(kappa^2/hy^2)*(y2(i,j,k)-y2(i,j-1,k))...
%                                     +(kappa^2/hz^2)*(y3(i,j,k)-y3(i,j,k-1));
%                                     %+imag(exp(-1i*y1(i,j,k))*conj(x(i,j,k))*x(i+1,j,k));
%                                     
%                     dPhidtY(i,j,k) = (kappa^2/hz^2)*(y2(i,j,k)+y2(i,j,k-1))...
%                                     +(kappa^2/hx^2)*(y2(i,j,k)+y2(i-1,j,k))...
%                                     +(kappa^2/hz^2)*(y3(i,j,k)-y3(i,j,k-1))...
%                                     +(kappa^2/hx^2)*(y1(i,j,k)-y1(i-1,j,k));
%                                     %+imag(exp(-1i*y2(i,j,k))*conj(x(i,j,k))*x(i,j+1,k));
%                                     
%                      dPhidtZ(i,j,k) = (kappa^2/hx^2)*(y3(i,j,k)+y3(i-1,j,k))...
%                                     +(kappa^2/hy^2)*(y3(i,j,k)+y3(i,j-1,k))...
%                                     +(kappa^2/hx^2)*(y1(i,j,k)-y1(i-1,j,k))...
%                                     +(kappa^2/hy^2)*(y2(i,j,k)-y2(i,j-1,k));
%                                     %+imag(exp(-1i*y3(i,j,k))*conj(x(i,j,k))*x(i,j,k+1));
