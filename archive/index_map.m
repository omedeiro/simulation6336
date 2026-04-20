function flat_index = index_map(i, j, k, p)
	if i < 1
		i = 1;
	elseif i > p.Nx
		i = p.Nx; 
	end

	if j < 1
		j = 1;
	elseif j > p.Ny
		j = p.Ny; 
	end

	if k < 1
		k = 1;
	elseif k > p.Nz
		k = p.z; 
	end

	flat_index = i + (j-1)*p.Ny + (k-1)*p.Nz^2;
end