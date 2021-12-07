function p = visualizeNetworkX2dInteract(t,X,p, app, event)

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

h1 = app.UIAxes;
% h1.Position = [900 200 550 800];

imagesc(app.UIAxes,'CData', abs(PSIT').^2./max(abs(PSIT').^2,1), 'HitTest','off')
% ax = gca;
set(app.UIAxes,'YDir','normal')
set(app.UIAxes, 'NextPlot','replacechildren')
colorbar(app.UIAxes)
caxis(app.UIAxes,[0 1])
colormap(app.UIAxes,'turbo')
axis(app.UIAxes, 'equal')
axis(app.UIAxes, 'tight')
% axis(app.UIAxes, 'off')
title(app.UIAxes,'|\psi|^2')
% xlabel(app.UIAxes,'x')
% ylabel(app.UIAxes,'y')
% zlabel(app.UIAxes,'z')


% app.UIFigureButtonDown(app);
pause(0.1)


Bxx  = reshape(Bx(ss2:se2,t-1),p.Nx-2,p.Ny-2);
Byy  = reshape(By(ss2:se2,t-1),p.Nx-2,p.Ny-2);
Bzz  = reshape(Bz(ss2:se2,t-1),p.Nx-2,p.Ny-2);

imagesc(app.UIAxes2,real(Bzz'))
set(app.UIAxes2,'YDir','normal')
colorbar(app.UIAxes2)

if max(max(abs(Bzz)))< 30e-3 && t>10/1e-1
    caxis(app.UIAxes2, [-30e-3 30e-3]);
else
    caxis(app.UIAxes2, [-1 1]);
end
colormap(app.UIAxes2,'turbo')
axis(app.UIAxes2, 'equal')
axis(app.UIAxes2, 'tight')
% axis(app.UIAxes2, 'off')
title(app.UIAxes2, 'B_z')




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


%%

end


