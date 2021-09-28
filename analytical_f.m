% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)

function F = analytical_f(x, hx, N)
    dPsidt = zeros(N, 1);
    dPhidt = zeros(N, 1);

    for i = 1:N
        if i==1
            dPsidt(i) =  (-2*x(i) + exp(-1i*x(N+i))*x(i+1))/hx^2 + (1-x(i)*conj(x(i)))*x(i);
            dPhidt(i) = imag(exp(-1i*x(N+i-1))*conj(x(i))*x(i+1));
        end
        if i > 1 && i < N
            dPsidt(i) = (exp(1i*x(N + i-1))*x(i-1) - 2*x(i) + exp(-1i*x(N+i))*x(i+1))/hx^2 + (1-x(i)*conj(x(i)))*x(i);
            dPhidt(i) = imag(exp(-1i*x(N+i-1))*conj(x(i))*x(i+1));
        end
        if i == N
            dPsidt(i) = (exp(1i*x(N + i-1))*x(i-1) - 2*x(i))/hx^2 + (1-x(i)*conj(x(i)))*x(i);
%             dPhidt(i) = imag(exp(-1i*x(N+i-1))*conj(x(i))*x(i+1));
        end
    end
    
    F = [dPsidt;dPhidt];
end