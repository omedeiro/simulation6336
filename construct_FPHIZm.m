function FPHIZ = construct_FPHIZm(x, y1, y2, y3, p)


mk = (p.Nx+1)*(p.Ny+1);
mj = (p.Nx+1);

m = p.M2;
M = p.m2;  
                
FPHIZ = sparse(M, 1, (p.kappa^2/p.hx^2).*(-y1(m+mk)+y1(m)+y1(m+mk-1)-y1(m-1))...
                    +(p.kappa^2/p.hy^2).*(-y2(m+mk)+y2(m)+y2(m+mk-mj)-y2(m-mj))...
                    +imag(exp(-1i*y3(m)).*conj(x(m)).*x(m+mk)) );          

end                            

