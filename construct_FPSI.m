function F = construct_FPSI(x)
    fun = @(a) (1-conj(a).*a).*a;
    F = fun(x);
end