function [L] = eval_L(x, a, p, u)
%EVAL_L equation 25
%   p = [h]
%   x = [x_1 x_2 x_3...]
N = length(x)
L = zeros(N, 1);
for i = 2:N-1
    L(i) = 1/p(1)^2 * (a(i-1)*x(i-1) - 2*x(i) + a(i+1)*x(i+1)); 
end
% Boundary conditions for x, y (z is different and needs to be periodic). Supercurrent across boundry == 0 (eq 35)
L(1) = L(2)*a(1);
L(N) = L(N-1)*conj(a(N-1));
end

