% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)

function F = analytical_f(x, hx, N)
    dPsidt = zeros(N, 1);
    dPhidt = zeros(N, 1);

%     for i = 1:N
%         if i==1
% %             dPsidt(i) = -(1-x(i)*conj(x(i)))*x(i);
%             BCPsi = x(i+1)*exp(-1i*x(N+i)); %psi0
% %             dPsidt(i) = (1-BCPsi*conj(BCPsi))*BCPsi;
%             dPsidt(i) = (exp(-1i*x(N+i))*x(i+1))/hx^2 + (1-BCPsi*conj(BCPsi))*BCPsi;
% 
%             dPhidt(i) = imag(exp(-1i*x(N+i-1))*conj(BCPsi)*x(i+1));
%         end
%         if i > 1 && i < N
%             dPsidt(i) = (exp(1i*x(N + i-1))*x(i-1) - 2*x(i) + exp(-1i*x(N+i))*x(i+1))/hx^2 + (1-x(i)*conj(x(i)))*x(i);
%             dPhidt(i) = imag(exp(-1i*x(N+i))*conj(x(i))*x(i+1));
%         end
%         if i == N
% %             dPsidt(i) = -(1-x(i)*conj(x(i)))*x(i);
%             BCPsi = x(i-1)*exp(1i*x(N+i-1)); %psiN
% %             dPsidt(i) = (1-BCPsi*conj(BCPsi))*BCPsi;
%             dPsidt(i) = (-2*BCPsi + exp(1i*x(N + i-1))*x(i-1) - 2*BCPsi)/hx^2 + (1-BCPsi*conj(BCPsi))*BCPsi;
% 
% %             dPhidt(i) = imag(exp(-1i*x(N+i-1))*conj(x(i))*x(i+1));
%         end
%     end
    for i = 1:N
        if i > 1 && i < N
            dPsidt(i) = (exp(1i*x(N + i-1))*x(i-1) - 2*x(i) + exp(-1i*x(N+i))*x(i+1))/hx^2 + (1-x(i)*conj(x(i)))*x(i);
            dPhidt(i) = imag(exp(-1i*x(N+i))*conj(x(i))*x(i+1));
        else
            %BCs
            i_m1 = N-1;
            i_p1 = 2;
            dPsidt(i) = (exp(1i*x(N + i_m1))*x(i_m1) - 2*x(i) + exp(-1i*x(N+i))*x(i_p1))/hx^2 + (1-x(i)*conj(x(i)))*x(i);
            dPhidt(i) = imag(exp(-1i*x(N+i))*conj(x(i))*x(i_p1));
        end
    end
    
    
    F = [dPsidt;dPhidt];
end