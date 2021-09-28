function Janalytic = eval_analytical_jac(x0, hx, N)
    psi = x0(1:N);
    phix = x0(N:2*N);
    % state = cat(1, psi, phix);

    Jpsidpsi = zeros(N, N);
    Jpsidphix = zeros(N, N);
    Jphixdpsi = zeros(N, N);
    Jphixdphix = zeros(N, N);

    for i = 1:N
       if i > 1
          Jpsidpsi(i, i-1) = exp(1i*phix(i-1))/hx^2;
       end
       Jpsidpsi(i, i) = -2/hx^2 + 1 - 2*psi(i)*conj(psi(i));
       if i < N
           Jpsidpsi(i, i+1) = exp(-1i*phix(i))/hx^2;
       end
    end

    for i = 1:N
       if i > 1
          Jpsidphix(i, i-1) = 1i*exp(1i*phix(i-1))*psi(i-1)/hx^2;
       end
       if i < N
          Jpsidphix(i, i+1) = 1i*exp(-1i*phix(i))*psi(i+1)/hx^2;
       end
    end

    for i = 1:N
       if i < N
          Jphixdpsi(i, i) = imag(exp(-1i*phix(i))*psi(i+1)); 
          Jphixdpsi(i, i+1) = imag(exp(-1i*phix(i)*conj(psi(i)))); 
       end
    end

    for i = 1:N
        if i < N
            Jphixdphix(i, i) = real(-exp(-1j*phix(i))*conj(psi(i))*psi(i+1)); 
        end
    end

    Janalytic = [Jpsidpsi Jpsidphix; Jphixdpsi Jphixdphix];
end



