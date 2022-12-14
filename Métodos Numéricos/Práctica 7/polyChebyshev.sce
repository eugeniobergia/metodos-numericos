function p = polyChebyshev(n)
    T = zeros(n+1,1);
    x = poly([0 1], "x", "c");
    T(1) = 1;
    T(2) = x;
    for i = 3:n+1
        T(i) = 2 * x * T(i-1) - T(i-2);
    end
    p = T(n+1)
endfunction

function root = rootChebyshev(n)
 root = zeros(n)
 for k = 1:n
     root(k) = cos(( (2*(k-1)+1) * %pi) /(2*n));
 end
 
endfunction
//
//root = rootChebyshev(3);
//p = polyChebyshev(3);
//disp(horner(p,root(1)))
