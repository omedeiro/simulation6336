function F = F_U(U, p, no_u)
    p.magBz = U;
    [u,p] = feval(p.eval_u, p.t, p.X0,p);
    F = feval(p.eval_f, p.X0, p, u);
end
