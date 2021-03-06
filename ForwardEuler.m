function X = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,timestep,visualize)
% uses Forward Euler to simulate states model dx/dt=f(x,p,u)
% from state x_start at time t_start
% until time t_stop, with time intervals timestep
% eval_f is a string including the name of the function that evaluates f(x,p,u)
% eval_u is a string including the name of the funciton that evaluates u(t)
% 
% X = ForwardEuler(eval_f,x_start,p,eval_u,t_start,t_stop,timestep,visualize)

% copyright Owen Medeiors, MIT 1969

X(:,1) = x_start;
t(1) = t_start;
% if visualize
%    visualizeResults(t,X,1,'.b');
% end
for n=1:ceil((t_stop-t_start)/timestep),
   %tic
   dt = min(timestep, (t_stop-t(n)));
   t(n+1)= t(n) + dt;
   [u, P] = feval(eval_u, t(n), X(:,n), p, n);
   p = P;
   p.Brealxt(n) = p.Brealx;
   p.Brealyt(n) = p.Brealy;
   p.Brealzt(n) = p.Brealz;
   f = feval(eval_f, X(:,n), p, u);
   X(:,n+1)= X(:,n) +  dt * f;
   if visualize
% %       visualizeResults(t,X,n+1,'.b');
        visualizeNetwork(n,X,p)
%         save('X.mat','X');
   end
   %t_step = toc
end
