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


tolrGCR   = 1e-4;  % convergence criteria on the GCR residual inside Newton
epsMF     = 1e-4;
errf	    = 1e-3;
errDeltax   = 1e-3;
relDeltax   = 1;     % note this is equivalent to NOT specifying it
MaxIter     = 20;

FiniteDifference=0;
visualizeGCR = 0;
eval_Jf = 0;

p.Breal=0;
for n=1:ceil((t_stop-t_start)/timestep) % TIME INTEGRATION LOOP
   dt = min(timestep, (t_stop-t(n)));
   t(n+1)= t(n) + dt;
   t(n)
   % Explicit solve for n+1 time
   [u,P] = feval(eval_u,t(n),X(:,n),p);
   p = P;
   p.Brealxt(n) = p.Brealx;
   p.Brealyt(n) = p.Brealy;
   p.Brealzt(n) = p.Brealz;

   f = feval(eval_f, X(:,n), p, u);
   Xpresent = X(:,n) +dt*f;
   gamma=X(:,n)+dt/2*f;

   % Newton loop to solve nonlinear equation
   [x,converged,errf_k,errDeltax_k,relDeltax_k,iterations] = NewtonGCRtrap(eval_f,Xpresent,p,u,errf,errDeltax,relDeltax,MaxIter,visualizeGCR,FiniteDifference,eval_Jf,tolrGCR,epsMF, gamma, dt);
   % Exit Newton loop

   X(:,n+1)= x ;
   
    y(n) = eval_y(X(:,n),p);
        figure(1000)
        plot(t(1:n), abs(y), 'o')

   if visualize
%       visualizeResults(t,X,n+1,'.b');
        visualizeNetwork(n,X,p)

   end
end

end