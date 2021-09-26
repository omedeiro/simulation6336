function [L] = eval_L(a, h)
%EVAL_L equation 25
%   p = [h]
%   x = [x_1 x_2 x_3...]
N = length(a);
L = sym(zeros(5));
for i = 2:N-1
    L(i, i-1:i+1) = [conj(a(i-1)), sym(-2), a(i+1)]; 
end
% Boundary conditions for x, y (z is different and needs to be periodic). Supercurrent across boundry == 0 (eq 35)
% L(1,:) = L(2)*a(1).*ones(N,1)
% L(N) = L(N-1)*conj(a(N-1));

L = L./(h^2);
end

