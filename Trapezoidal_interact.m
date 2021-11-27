function [X, p, f, errf_k,errDeltax_k] = Trapezoidal(eval_f,x_start,p,eval_u,t_start,t_stop,timestep,visualize)
% uses Forward Euler to simulate states model dx/dt=f(x,p,u)
% from state x_start at time t_start
% until time t_stop, with time intervals timestep
% eval_f is a string including the name of the function that evaluates f(x,p,u)
% eval_u is a string including the name of the funciton that evaluates u(t)
%
% X = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,timestep,visualize)

% copyright Luca Daniel, MIT 2018
X(:,1) = x_start;
t(1) = t_start;
% if visualize
% %    visualizeResults(t,X,1,'.b');
%     visualizeNetwork(t,X,p);
% end
global click

p.BXT = sparse(length(p.M2B),ceil((t_stop-t_start)/timestep));
p.BYT = sparse(length(p.M2B),ceil((t_stop-t_start)/timestep));
p.BZT = sparse(length(p.M2B),ceil((t_stop-t_start)/timestep));

tolrGCR   = 1e-3;  % convergence criteria on the GCR residual inside Newton
epsMF     = 1e-5;
errf	    = 1e-2;
errDeltax   = 1e-2;
relDeltax   = 1;     % note this is equivalent to NOT specifying it
MaxIter     = 20;

FiniteDifference=0;
visualizeGCR = 0;
eval_Jf = 0;
converged = 1;

p.Breal=0;

dt = timestep;
n = 1;
p.t(n) = 0;
% for n=1:ceil((t_stop-t_start)/timestep) % TIME INTEGRATION LOOP
while p.t(n) < t_stop % TIME INTEGRATION LOOP
    [u,P] = feval(eval_u,p.t(n),X(:,n),p,n);
    p = P;
    p.Brealxt(n) = p.Brealx;
    p.Brealyt(n) = p.Brealy;
    p.Brealzt(n) = p.Brealz;
    
    % Explicit solve for n+1 time
    f = feval(eval_f, X(:,n), p, u);
    Xpresent = X(:,n) +dt*f;
    gamma=X(:,n)+dt/2*f;
    
    % Newton loop to solve nonlinear equation
    [x,converged,errf_k,errDeltax_k,relDeltax_k,iterations] = NewtonGCRtrap(eval_f,Xpresent,p,u,errf,errDeltax,relDeltax,MaxIter,visualizeGCR,FiniteDifference,eval_Jf,tolrGCR,epsMF, gamma, dt);
    click = 0;
    % Exit Newton loop
    if converged
        X(:,n+1)= x ;
        p.t(n+1)= p.t(n) + dt;
%         p.t(n)
        
        n = n + 1;
%         dt = dt*2;
    else
%         p.t(n)
        dt = dt/10;
        disp("timestep reduced : " + dt)
    end
    
    
    
    if p.linearize == 1
        y(n) = eval_y(X(:,n),p);
        figure(1000)
        plot(t(1:n), abs(y), 'o')
    end
    
    if visualize
%       visualizeResults(t,X,n+1,'.b');
        p = visualizeNetworkX2dInteract(n,X,p);
        
    end
end

end