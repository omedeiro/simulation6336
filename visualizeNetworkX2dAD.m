function h1 = visualizeNetworkX2d(X,p)


%     [Bx,By,Bz] = eval_Bfield(X,p);
    Bx = p.BXT;
    By = p.BYT;
    Bz = p.BZT;
   
    [yy, xx] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-1), p.hy:p.hy:p.hy*(p.Ny-1));

    [yy2, xx2] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-2), p.hy:p.hy:p.hy*(p.Ny-2));
    xx = cube2column(xx);
    yy = cube2column(yy);
    xx2 = cube2column(xx2);
    yy2 = cube2column(yy2);
    mk = (p.Nx-1)*(p.Ny-1); 
    mk2 = (p.Nx-2)*(p.Ny-2); 

    n = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
    n2 = (p.Nx-2)*(p.Ny-2)*(p.Nz-2);
    
    
    % Slice
    if p.slice > 0
        ss = mk * (p.slice - 1) + 1;
        se   = mk * (p.slice);
        ss2 = mk2 * (p.slice - 1) + 1;
        se2   = mk2 * (p.slice);
    else
        ss = 1;
        se = n;
        ss2 = 1;
        se2 = n;
    end
        
    
for t = floor(linspace(1,size(X,2)-1,p.frames))
    
    PSIT = X(1:n,t);
    PHITX = X(n+1:2*n,t);
    PHITY = X(2*n+1:3*n,t);
    PHITZ = X(3*n+1:4*n,t);

    PSIT = PSIT(ss:se);
    PHITX = PHITX(ss:se);
    PHITY = PHITY(ss:se);
    PHITZ = PHITZ(ss:se);
    xx = cube2column(xx);
    yy = cube2column(yy);
    PSIT = reshape(PSIT,p.Nx-1,p.Ny-1);
    PHITX = reshape(PHITX,p.Nx-1,p.Ny-1);
    PHITY = reshape(PHITY,p.Nx-1,p.Ny-1);
    PHITZ = reshape(PHITZ,p.Nx-1,p.Ny-1);

    h1 = figure(1);
%     h1.Position = [900 200 550 800];
    
    subplot(5,3,1:9)
    imagesc(xx,yy,abs(PSIT).^2./max(abs(PSIT).^2,1))
    colorbar
    caxis([0 1])
    axis equal
    axis tight
    axis off
    title('|\psi|^2')
    xlabel('x')
    ylabel('y')
    zlabel('z')


    subplot(5,3,10)
    imagesc(xx,yy,real(PHITX))
    colorbar
    title('\phi_x')
    axis equal
    axis tight
    axis off
    xlabel('x')
    ylabel('y')
    zlabel('z')

    subplot(5,3,11)
    imagesc(xx,yy,real(PHITY))
    colorbar
    title('\phi_y')
    axis equal
    axis tight
    axis off
    xlabel('x')
    ylabel('y')
    zlabel('z')

    subplot(5,3,12)
    imagesc(xx,yy,real(PHITZ))
    colorbar
    title('\phi_z')
    axis equal
    axis tight
    axis off
    xlabel('x')
    ylabel('y')
    zlabel('z')






    Bxx  = reshape(Bx(ss2:se2,t),p.Nx-2,p.Ny-2);
    Byy  = reshape(By(ss2:se2,t),p.Nx-2,p.Ny-2);
    Bzz  = reshape(Bz(ss2:se2,t),p.Nx-2,p.Ny-2);

    subplot(5,3,13)
    imagesc(xx2,yy2,real(Bxx))
    colorbar
    axis equal
    axis tight
    axis off
    title('B_x')
    xlabel('x')
    ylabel('y')
    zlabel('z')


    subplot(5,3,14)
    imagesc(xx2,yy2,real(Byy))
    colorbar
    axis equal
    axis tight
    axis off
    title('B_y')
    xlabel('x')
    ylabel('y')
    zlabel('z')

    subplot(5,3,15)
    imagesc(xx2,yy2,real(Bzz))
    colorbar
%     hold on
%     contour(real(Bzz), [-0.02 0 0.02],'--', 'LineWidth', 2)
    axis equal
    axis tight
    axis off
    title('B_z')
    xlabel('x')
    ylabel('y')
    zlabel('z')

    sgtitle({"applied Bx = "+p.Brealxt(t)+", applied By = "+p.Brealyt(t)+", applied Bz = "+p.Brealzt(t),"t = "+ p.t(t),})

%     drawnow
%     if p.visualizeSave == 1
%         if t ==1
%             addpath(pwd+"/gifs")
%             filename = char(pwd+"/gifs/trapGif"+datestr(now,30)+".gif");
%             gif(filename);
%             save(filename(1:end-4),'p')
%         else
%             gif
%         end
%     end

end

%%

end


