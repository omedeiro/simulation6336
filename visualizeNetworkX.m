function visualizeNetworkX(X,p)

[Bx,By,Bz] = eval_Bfield(X,p);

[yy, xx, zz] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-1), p.hy:p.hy:p.hy*(p.Ny-1), p.hz:p.hz:p.hz*(p.Nz-1));

[yy2, xx2, zz2] = meshgrid(p.hx:p.hx:p.hx*(p.Nx-2), p.hy:p.hy:p.hy*(p.Ny-2), p.hz:p.hz:p.hz*(p.Nz-2));
xx = cube2column(xx);
yy = cube2column(yy);
zz = cube2column(zz);
xx2 = cube2column(xx2);
yy2 = cube2column(yy2);
zz2 = cube2column(zz2);
n = (p.Nx-1)*(p.Ny-1)*(p.Nz-1);
n2 = (p.Nx-2)*(p.Ny-2)*(p.Nz-2);
for t = floor(linspace(1,size(X,2),50))
omm=1;
if omm ==1
    figure(3)
    PSIT = X(1:n,t);
    PHITX = X(n+1:2*n,t);
    PHITY = X(2*n+1:3*n,t);
    PHITZ = X(3*n+1:4*n,t);

    subplot(2,2,1)
    scatter3(xx,yy,zz,36,(abs(PSIT).^2)./max(max((abs(PSIT).^2))), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
%     caxis([0 1])
    axis equal
    title({'|\psi|^2',''})
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view(0,90)

    subplot(2,2,2)
    scatter3(xx2,yy2,zz2,36,real(Bx(:,t)), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
%     caxis([-1 1])
    axis equal
    title({'B_x',''})
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view(0,90)

        
    subplot(2,2,3)
    scatter3(xx2,yy2,zz2,36,real(By(:,t)), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
%     caxis([-1 1])
    axis equal
    title({'B_y',''})
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view(0,90)

    subplot(2,2,4)
    scatter3(yy2,xx2,zz2,36,real(Bz(:,t)), 'filled', 'MarkerFaceAlpha', 0.5)

%         scatter3(yy2,xx2,zz2,36,real(Bz(:,t))./max(max(real(Bz(:,t))),1), 'filled', 'MarkerFaceAlpha', 0.5)
    colorbar
%     caxis([-1 1])
    axis equal
    title({'B_z',''})
    xlabel('x')
    ylabel('y')
    zlabel('z')
    view(0,90)
    sgtitle("t = "+ t*p.timestep + ", applied Bz = "+p.Brealt(t))
    pause(0.005)
    drawnow
    if p.visualizeSave == 1
        if t ==1
            addpath(pwd+"/gifs")
            filename = char(pwd+"/gifs/trapGif"+datestr(now,30)+".gif");
            gif(filename);
            save(filename(1:end-4),'p')
        else
            gif
        end
    end
    % Capture the plot as an image 
%     frame = getframe(gcf); 
%     im = frame2im(frame); 
%     [imind,cm] = rgb2ind(im,256); 
%     % Write to the GIF File 
%     if t == 1 
%       imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
%     else 
%       imwrite(imind,cm,filename,'gif','WriteMode','append'); 
%     end 
end
    %%

    Bavg_x = real(mean(Bx(:,t),'all', 'omitnan'));
    Bavg_y = real(mean(By(:,t),'all', 'omitnan'));
    Bavg_z = real(mean(Bz(:,t),'all', 'omitnan'));
    
    om = 0;
    if om==1
    figure(3)
    q = quiver3(1,1,1,Bavg_x,Bavg_y,Bavg_z);
    xlim([0 2])
    ylim([0 2])
    end
    
    %%
    ommm = 0;
    if ommm==0
    figure(4)
        PSIT = X(1:n,t);
        PHITX = X(n+1:2*n,t);
        PHITY = X(2*n+1:3*n,t);
        PHITZ = X(3*n+1:4*n,t);

        subplot(2,2,1)
        scatter3(xx,yy,zz,36,abs(PSIT).^2./max(max(abs(PSIT).^2),1), 'filled', 'MarkerFaceAlpha', 0.5)
        colorbar
        caxis([0 1])
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')

        subplot(2,2,2)
        scatter3(xx,yy,zz,36,abs(PHITX).^2./max(max(abs(PHITX).^2),1), 'filled', 'MarkerFaceAlpha', 0.5)
        colorbar
        caxis([0 1])
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')

        subplot(2,2,3)
        scatter3(xx,yy,zz,36,abs(PHITY).^2./max(max(abs(PHITY).^2),1), 'filled', 'MarkerFaceAlpha', 0.5)
        colorbar
        caxis([0 1])
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')

        subplot(2,2,4)
        scatter3(xx,yy,zz,36,abs(PHITZ).^2./max(max(abs(PHITZ).^2),1), 'filled', 'MarkerFaceAlpha', 0.5)
        colorbar
        caxis([0 1])
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')

        pause(0.005)    
        sgtitle("t = "+ t*p.timestep + ", applied Bz = "+p.Brealt(t))
        drawnow
    end

end
end