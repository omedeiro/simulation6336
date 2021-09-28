
psi = [1; 0; 2];
phix = [0; 1; 2];

hx = 1;
% state = cat(1, psi, phix);

Jpsidpsi = zeros(N, N);
Jpsidphix = zeros(N, N);
Jphixdpsi = zeros(N, N);
Jphixdphix = zeros(N, N);

i_m1 = N-1;
i_p1 = 2;

% for i = 1:N
%    if i > 1
%       Jpsidpsi(i, i-1) = exp(1i*phix(i-1))/hx^2;
%    end
%    Jpsidpsi(i, i) = -2/hx^2 + psi(i)^2 + 2*psi(i)*conj(psi(i));
%    if i < N
%        Jpsidpsi(i, i+1) = exp(-1i*phix(i))/hx^2;
%    end
% end

for i = 1:N
   if i > 1 && i < N
      Jpsidpsi(i, i-1) = exp(1i*phix(i-1))/hx^2;
      Jpsidpsi(i, i+1) = exp(-1i*phix(i))/hx^2;
   else
      %BCs
      Jpsidpsi(i, i_m1) = exp(1i*phix(i_m1))/hx^2;
      Jpsidpsi(i, i_p1) = exp(-1i*phix(i))/hx^2;
   end
   Jpsidpsi(i, i) = -2/hx^2 +1 -2*psi(i)*conj(psi(i));
end


for i = 1:N
   if i > 1 && i < N
      Jpsidphix(i, i-1) = 1i*exp(1i*phix(i-1))*psi(i-1)/hx^2;
      Jpsidphix(i, i+1) = -1i*exp(-1i*phix(i))*psi(i+1)/hx^2;
   else
      %BCs
      Jpsidphix(i, i_m1) = 1i*exp(1i*phix(i_m1))*psi(i_m1)/hx^2;
      Jpsidphix(i, i_p1) = -1i*exp(-1i*phix(i))*psi(i_p1)/hx^2;
   end
end


% for i = 1:N
%    if i < N
%       Jphixdpsi(i, i) = imag(exp(1i*phix(i))*psi(i+1)); 
%       Jphixdpsi(i, i+1) = imag(exp(1i*phix(i)*conj(psi(i)))); 
%    end
% end

for i = 1:N
    if i ~= N
        Jphixdpsi(i, i) = imag(exp(1i*phix(i))*psi(i+1)); 
        Jphixdpsi(i, i+1) = imag(exp(1i*phix(i)*conj(psi(i)))); 
    else
        % BCs
        Jphixdpsi(i, i) = imag(exp(1i*phix(i))*psi(i_p1)); 
        Jphixdpsi(i, i_p1) = imag(exp(1i*phix(i)*conj(psi(i)))); 
    end
end


for i = 1:N
    if i ~= N
        Jphixdphix(i, i) = real(-exp(-1i*phix(i))*conj(psi(i))*psi(i+1)); 
    else
        % BCs
        Jphixdphix(i, i) = real(-exp(-1i*phix(i))*conj(psi(i))*psi(i_p1));
    end
end

Janalytic = [Jpsidpsi Jpsidphix; Jphixdpsi Jphixdphix];


