function F_PSI = construct_FPSI(x)
    fun = @(a) (1-conj(a).*a).*a;
    F_PSI = fun(x);
end