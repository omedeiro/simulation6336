function visualizeNetworkFace(X,p)

figure(10)
for n=1:size(X,2)
   if n==1
        psi = X(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),n+1);
        psi2 = column2cube(abs(psi).^2, (p.Nx-1), (p.Ny-1), (p.Nz-1));
        psi_surf = psi2(:,:,1);
        s = surf(psi_surf)
        view(0,90)
        colorbar
        caxis([0 max(max(psi_surf))])

   else
       psi = X(1:(p.Nx-1)*(p.Ny-1)*(p.Nz-1),n);
       psi2 = column2cube(abs(psi).^2, (p.Nx-1), (p.Ny-1), (p.Nz-1));
       s.ZData = psi2(:,:,1);
       caxis([0 max(max(psi_surf))])

   end
   pause(0.05)
end

end