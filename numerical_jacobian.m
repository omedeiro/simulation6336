N = 3;
psi = [1; 1; 1];
phix = [1; 1; 1];

hx = 1;

state_0 = cat(1, psi, phix);

%define F

J = eval_num_jac(state_0, F);


function J = eval_num_jac(x0, F)
    eps_Im = 0.1;
    eps_Re = 0.1;
    S_Im = 2;
    S_Re = 2;
    J = zeros(size(F(x0),1), size(x0,1));
    Jp = ones(size(F(x0),1), size(x0,1));
    err = 1e-16;
    while any(abs((J - Jp)) > err, 'all') 
        Jp = J;
        for k = 1 : size(J,2)
            dx = zeros(size(x0,1), 1);
            dx(k) = eps_Re + eps_Imi;
            J(:, k) = (F(x0 + dx) - F(x0))/eps;
        end
        eps_Re = eps_Re/S_Re;
        eps_Im = eps_Im/S_Im;
    end
end