% x is state variable fxn x = [psi1, psi2, ...psiN, phix1, phix2 ... phixN]
% --> phix1 = x(N+1)
% 
% x(i) = x(i+1)*exp(-1i*x(N+i)); %psi1
% x(N) = x(i-1)*exp(1i*x(N+i-1)); %psiN
function Janalytic = eval_analytical_jac_xyz(psi, phix, phiy, phiz, hx, hy, hz, kappa, N)

    Jpsidpsi = zeros(N^3, N^3);
    Jpsidphix = zeros(N^3, N^3);
    Jpsidphiy = zeros(N^3, N^3);
    Jpsidphiz = zeros(N^3, N^3);
    
    Jphixdpsi = zeros(N^3, N^3);
    Jphixdphix = zeros(N^3, N^3);
    Jphixdphiy = zeros(N^3, N^3);
    Jphixdphiz = zeros(N^3, N^3);
    
    Jphiydpsi = zeros(N^3, N^3);
    Jphiydphix = zeros(N^3, N^3);
    Jphiydphiy = zeros(N^3, N^3);
    Jphiydphiz = zeros(N^3, N^3);
    
    Jphizdpsi = zeros(N^3, N^3);
    Jphizdphix = zeros(N^3, N^3);
    Jphizdphiy = zeros(N^3, N^3);
    Jphizdphiz = zeros(N^3, N^3);
    
    for k = 1:N
    	for j = 1:N
    		for i = 1:N
    			current = index_map(i, j, k, N);
    			
    			% derivatives of psi wrt psi
    			Jpsidpsi(current, current) = 2/hx^2 + 2/hy^2 + 2/hz^2 - 2*psi(current)*conj(psi(current));
    			
    			Jpsidpsi(current, index_map(i-1, j, k, N)) = exp(1i*phix(index_map(i-1, j, k, N)))/hx^2;
    			Jpsidpsi(current, index_map(i+1, j, k, N)) = exp(-1i*phix(current))/hx^2;
    			
    			Jpsidpsi(current, index_map(i, j-1, k, N)) = exp(1i*phiy(index_map(i, j-1, k, N)))/hy^2;
    			Jpsidpsi(current, index_map(i, j+1, k, N)) = exp(-1i*phiy(current))/hy^2;
    			
    			Jpsidpsi(current, index_map(i, j, k-1, N)) = exp(1i*phiz(index_map(i, j, k-1, N)))/hz^2;
    			Jpsidpsi(current, index_map(i, j, k+1, N)) = exp(-1i*phiz(current))/hz^2;
    			
    			%derivatives of psi wrt phi
    			Jpsidphix(current, current) = -1i*exp(-1i*phix(current))*psi(index_map(i+1, j, k, N))/hx^2;
    			Jpsidphix(current, index_map(i-1, j, k, N)) = 1i*exp(1i*phix(index_map(i-1, j, k, N)))*psi(index_map(i-1, j, k, N))/hx^2;
    			
    			Jpsidphiy(current, current) = -1i*exp(-1i*phiy(current))*psi(index_map(i, j+1, k, N))/hy^2;
    			Jpsidphiy(current, index_map(i, j-1, k, N)) = 1i*exp(1i*phiy(index_map(i, j-1, k, N)))*psi(index_map(i, j-1, k, N))/hy^2;
    			
    			Jpsidphiz(current, current) = -1i*exp(-1i*phiz(current))*psi(index_map(i, j, k+1, N))/hz^2;
    			Jpsidphiz(current, index_map(i, j, k-1, N)) = 1i*exp(1i*phiz(index_map(i, j, k-1, N)))*psi(index_map(i, j, k-1, N))/hz^2;
    			
    			% derivatives of phix
    			[diagonal, right] = dphidpsi(phix, psi, hx, i, j, k, N);
    			Jphixdpsi(current, current) = diagonal;
    			Jphixdpsi(current, index_map(i+1, j, k, N)) = right;
    			
    			[diagonal, offy, offz] = dphidself(kappa, hy, hz);
    			Jphixdphix(current, current) = diagonal;
    			Jphixdphix(current, index_map(i, j+1, k, N)) = offy;
    			Jphixdphix(current, index_map(i, j-1, k, N)) = offy;
    			Jphixdphix(current, index_map(i, j, k+1, N)) = offz;
    			Jphixdphix(current, index_map(i, j, k-1, N)) = offz;
    			
    			[positive, negative] = dphindphim(kappa, hy);
    			Jphixdphiy(current, current) = positive;
    			Jphixdphiy(current, index_map(i+1, j, k, N)) = negative;
    			Jphixdphiy(current, index_map(i+1, j-1, k, N)) = positive;
    			Jphixdphiy(current, index_map(i, j-1, k, N)) = negative;
    			
    			[positive, negative] = dphindphim(kappa, hz);
    			Jphixdphiz(current, current) = positive;
    			Jphixdphiz(current, index_map(i+1, j, k, N)) = negative;
    			Jphixdphiz(current, index_map(i+1, j, k-1, N)) = positive;
    			Jphixdphiz(current, index_map(i, j, k-1, N)) = negative;
    			
    			% derivatives of phiy
    			[diagonal, right] = dphidpsi(phiy, psi, hy, i, j, k, N);
    			Jphiydpsi(current, current) = diagonal;
    			Jphiydpsi(current, index_map(i, j+1, k, N)) = right;
    			
    			[diagonal, offz, offx] = dphidself(kappa, hz, hx);
    			Jphiydphiy(current, current) = diagonal;
    			Jphiydphiy(current, index_map(i, j, k+1, N)) = offz;
    			Jphiydphiy(current, index_map(i, j, k-1, N)) = offz;
    			Jphiydphiy(current, index_map(i-1, j, k, N)) = offx;
    			Jphiydphiy(current, index_map(i+1, j, k, N)) = offx;
    			
    			[positive, negative] = dphindphim(kappa, hx);
    			Jphiydphix(current, current) = positive;
    			Jphiydphix(current, index_map(i, j+1, k, N)) = negative;
    			Jphiydphix(current, index_map(i-1, j+1, k, N)) = positive;
    			Jphiydphix(current, index_map(i-1, j, k, N)) = negative;
    			
    			[positive, negative] = dphindphim(kappa, hz);
    			Jphiydphiz(current, current) = positive;
    			Jphiydphiz(current, index_map(i, j+1, k, N)) = negative;
    			Jphiydphiz(current, index_map(i, j+1, k-1, N)) = positive;
    			Jphiydphiz(current, index_map(i, j, k-1, N)) = negative;
    			
    			% derivatives of phiz
    			[diagonal, right] = dphidpsi(phiz, psi, hz, i, j, k, N);
    			Jphizdpsi(current, current) = diagonal;
    			Jphizdpsi(current, index_map(i, j, k+1, N)) = right;
    			
    			[diagonal, offx, offy] = dphidself(kappa, hx, hy);
    			Jphizdphiz(current, current) = diagonal;
    			Jphizdphiz(current, index_map(i+1, j, k, N)) = offx;
    			Jphizdphiz(current, index_map(i-1, j, k, N)) = offx;
    			Jphizdphiz(current, index_map(i, j+1, k, N)) = offy;
    			Jphizdphiz(current, index_map(i, j-1, k, N)) = offy;
    			
    			[positive, negative] = dphindphim(kappa, hx);
    			Jphizdphix(current, current) = positive;
    			Jphizdphix(current, index_map(i, j, k+1, N)) = negative;
    			Jphizdphix(current, index_map(i-1, j, k+1, N)) = positive;
    			Jphizdphix(current, index_map(i-1, j, k, N)) = negative;
    			
    			[positive, negative] = dphindphim(kappa, hy);
    			Jphizdphiy(current, current) = positive;
    			Jphizdphiy(current, index_map(i, j, k+1, N)) = negative;
    			Jphizdphiy(current, index_map(i, j-1, k+1, N)) = positive;
    			Jphizdphiy(current, index_map(i, j-1, k, N)) = negative;
    		end 
    	end
    end
% 	spy(Jpsidpsi)
    Janalytic = [Jpsidpsi Jpsidphix Jpsidphiy Jpsidphiz; Jphixdpsi Jphixdphix Jphixdphiy Jphixdphiz; Jphiydpsi Jphiydphix Jphiydphiy Jphiydphiz; Jphizdpsi Jphizdphix Jphizdphiy Jphizdphiz];
end


function flat_index = index_map(i, j, k, N)
	if i < 1
		i = 1;
	elseif i > N
		i = N; 
	end

	if j < 1
		j = 1;
	elseif j > N
		j = N; 
	end

	if k < 1
		k = 1;
	elseif k > N
		k = N; 
	end

	flat_index = i + (j-1)*N + (k-1)*N^2;
end

function [diagonal, right] = dphidpsi(phi, psi, h, i, j, k, N)
	current = index_map(i, j, k, N);
	diagonal = imag(exp(-1i*phi(current))*psi(index_map(i+1, j, k, N)));
	right = imag(exp(-1i*phi(current)*conj(psi(current))));
end

function [diagonal, offdiagy, offdiagz] = dphidself(kappa, hy, hz)
	diagonal = -2*kappa^2/hy^2 - 2*kappa^2/hz^2;
	offdiagy = kappa^2/hy^2;
	offdiagz = kappa^2/hz^2;
end

function [positive, negative] = dphindphim(kappa, hm)
	positive = kappa^2/hm^2;
	negative = -kappa^2/hm^2;
end