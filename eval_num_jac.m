function J = eval_num_jac(x0, F)
    eps_Im = .01;
    eps_Re = .01;
    S_Im = 2;
    S_Re = 2;
    J = zeros(size(F(x0),1), size(x0,1));
    Jp = ones(size(F(x0),1), size(x0,1));
    err = 1e-8;
    while any(abs((J - Jp)) > err, 'all') 
        Jp = J;
        for k = 1 : size(J,2) % loop columns
            dx = zeros(size(x0,1), 1);
            dx(k) = eps_Re + 1i*eps_Im;
            J(:, k) = (F(x0 + dx) - F(x0))/dx(k);
        end
        eps_Re = eps_Re/S_Re;
        eps_Im = eps_Im/S_Im;
    end
end
