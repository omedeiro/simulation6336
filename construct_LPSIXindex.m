function LPSI = construct_LPSIXindex(y,p)
zzz = 0;
LPSI = spalloc((p.Nx+1)*(p.Ny+1)*(p.Nz+1), (p.Nx+1)*(p.Ny+1)*(p.Nz+1), (p.Nx+1)*(p.Ny+1)*(p.Nz+1));
    i = 1;    
    j = 1;
    k = 1;
    
    for z = 1:length(y)

        m = (p.Nx+1)*(j-1)+(p.Nx+1)*(p.Ny+1)*(k-1);
            
        if i~=1 && j~=1 && k~=1 
            LPSI(z+m,z-1+m) = exp(-1i*y(z));
            LPSI(z+m,z+m) = -2;
            LPSI(z+m,z+1+m) = exp(1i*y(z));
            
        end
        
        if i == p.Nx && j ~= p.Ny
            i = 1;
            j= j+1;

        elseif i == p.Nx && j == p.Ny
            i = 1;
            j = 1;
            k = k+1;
        else
            i=i+1;
        end
        zzz = zzz+1;
        disp(zzz)
    end
end