function [F,p] = eval_linf(X, p, U)
    
    F = p.A*X + p.B*[1; U];
end