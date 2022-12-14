function i = trapecio(f, a, b)
    h = b - a;
    i = (f(a) + f(b))*h/2;
endfunction
/*
--> trapecio(exp, 0, 3)
 ans  =

   31.628305
*/

function i = trapecio_compuesto(f, a, b, n)
    h = (b - a)/n;
    i = f(a)/2 + f(b)/2;
    for j = 1:(n - 1)
        i = i + f(a + j*h);
    end
    i = i*h;
endfunction
/*
--> trapecio_compuesto(exp, 0, 3, 10)
 ans  =

   19.228464
*/

function i = simpson(f, a, b)
    x1 = (a + b)/2;
    h = b - x1;
    i = (f(a) + 4*f(x1) + f(b))*h/3;
endfunction
/*
--> simpson(exp, 0, 3)
 ans  =

   19.506147
*/

function i = simpson_compuesto(f, a, b, n)
    _n = n/2;
    if int(_n) <> _n then
        error("simpson_compuesto - n debe ser par.");
        abort();
    end
    h = (b - a)/n;
    i = f(a) + f(b);
    for j = 1:(n - 1)
        _j = j/2;
        if int(_j) == _j then
            i = i + 2*f(a + j*h);
        else
            i = i + 4*f(a + j*h);
        end
    end
    i = i*h/3;
endfunction
/*
--> simpson_compuesto(exp, 0, 3, 10)
 ans  =

   19.086387
*/

function e = err_trapecio(a, b, n, cota)
    h = (b - a)/n;
    e = -(h**2)*(b - a)*cota/12;
endfunction
/*
--> err_trapecio(0, 3, 10, exp(3))
 ans  =

  -0.4519246
*/

function e = err_simpson(a, b, n, cota)
    h = (b - a)/n;
    e = -(h**4)*(b - a)*cota/180;
endfunction
/*
--> err_simpson(0, 3, 10, exp(3))
 ans  =

  -0.0027115
*/
