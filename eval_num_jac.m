function J = eval_num_jac(X0, F, err)
    eps_Im = .01;
    eps_Re = .01;
    S_Im = 2;
    S_Re = 2;
    J = zeros(size(F(X0),1), numel(X0));
    Jp = ones(size(F(X0),1), numel(X0));
    while any(abs((J - Jp)) > err, 'all') 
        Jp = J;
        for k = 1 : size(J,2) % loop columns
            dx = zeros(size(X0,1), 1);
            dx(k) = eps_Re + 1i*eps_Im;
            J(:, k) = (F(X0 + dx) - F(X0))/dx(k);
        end
        eps_Re = eps_Re/S_Re;
        eps_Im = eps_Im/S_Im;
        disp(max(max(abs((J - Jp)))))
    end
end

