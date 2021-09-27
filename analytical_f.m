% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)

function dpsidt = evolve_psi(x, hx, N)
    dpsidt = zeros(1, N);
    
    for i = 1:N
        dpsidt(1, i) = (exp(1i*x(N + i-1))*x(i-1) - 2*x(i) + exp(-1i*x(N+i))*x(i+1))/hx^2 + (1-x(i)*conj(x(i)))*x(i);
    end
end