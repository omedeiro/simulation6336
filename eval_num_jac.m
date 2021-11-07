function J = eval_num_jac(X0, eval_f, p, eval_u, t, err)
    eps_Im = .01;
    eps_Re = .01;
    S_Im = 2;
    S_Re = 2;
    u = feval(eval_u,t);
   
    J = sparse(size(feval(eval_f, X0, p, u),1), numel(X0));
    Jp = sparse(size(feval(eval_f, X0, p, u),1), numel(X0),1);
    while any(abs((J - Jp)) > err, 'all') 
        Jp = J;
        for k = 1 : size(J,2) % loop columns
            dx = sparse(size(X0,1),1);
            dx(k) = eps_Re + 1i*eps_Im;
            J(:, k) = (feval(eval_f, X0+dx, p, u) - feval(eval_f, X0, p, u))/dx(k);
        end
        eps_Re = eps_Re/S_Re;
        eps_Im = eps_Im/S_Im;
        disp(max(max(abs((J - Jp)))))
    end
end

