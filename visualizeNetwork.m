function visualizeNetwork(X,p)

[Bx,By,Bz] = eval_Bfield(X,p);

[xx, yy, zz] = meshgrid(0:p.hx:p.hx*(p.Nx-2), 0:p.hy:p.hy*(p.Ny-2), 0:p.hz:p.hz*(p.Nz-2));

[xx2, yy2, zz2] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-2), p.hy:p.hy:p.hy*(p.Ny-2), p.hz:p.hz:p.hz*(p.Nz-2));
xx = cube2column(xx);
yy = cube2column(yy);
zz = cube2column(zz);
xx2 = cube2column(xx2);
yy2 = cube2column(yy2);
zz2 = cube2column(zz2);
n = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
n2 = (p.Nx-2)*(p.Ny-2)*(p.Nz-2);
for t = 1:size(X, 2)
    figure(2)
    PSIT = X(1:n,t);
    PHITX = X(n+1:2*n,t);
    PHITY = X(2*n+1:3*n,t);
    PHITZ = X(3*n+1:4*n,t);

    subplot(2,2,1)
    scatter3(xx,yy,zz,36,abs(PSIT), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    caxis([0 max(max(abs(PSIT)))]);
    axis equal
    title('\psi')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,2)
    scatter3(xx2,yy2,zz2,36,abs(Bx(:,t)), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    caxis([0 max(max(Bx))]);
    axis equal
    title('Bx')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,3)
    scatter3(xx2,yy2,zz2,36,abs(By(:,t)), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    caxis([0 max(max(By))]);
    axis equal
    title('By')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    subplot(2,2,4)
    scatter3(xx2,yy2,zz2,36,abs(Bz(:,t)), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
    caxis([0 max(max(Bz))]);
    axis equal
    title('Bz')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    
    pause(0.5)    
    sgtitle(t)
    drawnow
    
%     figure(3)
%     PSIT = X(1:n,t);
%     PHITX = X(n+1:2*n,t);
%     PHITY = X(2*n+1:3*n,t);
%     PHITZ = X(3*n+1:4*n,t);
% 
%     subplot(2,2,1)
%     scatter3(xx,yy,zz,36,imag(PSIT), 'filled', 'MarkerFaceAlpha', 0.5)
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     
%     subplot(2,2,2)
%     scatter3(xx,yy,zz,36,imag(PHITX), 'filled', 'MarkerFaceAlpha', 0.5)
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     
%     subplot(2,2,3)
%     scatter3(xx,yy,zz,36,real(PHITY), 'filled', 'MarkerFaceAlpha', 0.5)
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     
%     subplot(2,2,4)
%     scatter3(xx,yy,zz,36,real(PHITZ), 'filled', 'MarkerFaceAlpha', 0.5)
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     
%     pause(0.05)    
%     sgtitle(t)
%     drawnow
%     
%     figure(4)
%     subplot(1,3,1)
%     scatter3(xx2,yy2,zz2,36,real(Bx(:,t)),'filled', 'MarkerFaceAlpha',0.5)
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     
%     subplot(1,3,2)
%     scatter3(xx2,yy2,zz2,36,real(By(:,t)),'filled', 'MarkerFaceAlpha',0.5)
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     
%     subplot(1,3,3)
%     scatter3(xx2,yy2,zz2,36,real(Bz(:,t)),'filled', 'MarkerFaceAlpha',0.5)
%     colorbar
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
end

end