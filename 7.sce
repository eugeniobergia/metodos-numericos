function l = L_i(p, x, i)
    l = 1;
    n = length(x);
    
    for j = 1:n
        if j <> i then
            l = l * (p - x(j)) / (x(i) - x(j));
        end
    end
endfunction
/*
--> x = [0 0.2 0.4 0.6]';
--> L_i(poly([0], "x"), x, 3)
 ans  =

            2       3
  -7.5x +50x  -62.5x
*/

function l = lagrange(p, x, y)
    l = 0;
    n = length(x);
    
    for i = 1:n
        l = l + y(i) * L_i(p, x, i);
    end
endfunction
/*
--> x = [0 0.2 0.4 0.6]';
--> y = [1.0 1.2214 1.4918 1.8221]';
--> lagrange(poly([0], "x"), x, y)
 ans  =

                          2            3
   1 +1.0026667x +0.47625x  +0.2270833x
*/

function d = dd(x, y)
    n = length(x);
    
    if n == 2 then
        d = (y(n) - y(1)) / (x(n) - x(1));
    else
        d = (dd(x(2:n), y(2:n)) - dd(x(1:n-1), y(1:n-1))) / (x(n) - x(1));
    end
endfunction
/*
--> x = [0 0.2 0.4 0.6]';
--> y = [1.0 1.2214 1.4918 1.8221]';
--> dd(x, y)
 ans  =

   0.2270833
*/

function d = newton_dd(p, x, y)
    d = 0;
    n = length(x);
    
    for i = n:-1:2
        d = (d + dd(x(1:i), y(1:i))) * (p - x(i-1));
    end
    
    d = d + y(1);
endfunction
/*
--> x = [0 0.2 0.4 0.6]';
--> y = [1.0 1.2214 1.4918 1.8221]';
--> newton_dd(poly([0], "x"), x, y)
 ans  =

                          2            3
   1 +1.0026667x +0.47625x  +0.2270833x
*/

function [x, a, it] = gausselimPP(A, b)
    [nA, mA] = size(A)
    [nb, mb] = size(b)
    
    if nA <> mA then
        error('gausselimPP - La matriz A debe ser cuadrada');
        abort;
    elseif mA <> nb then
        error('gausselimPP - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    it = 0;
    a = [A b]; // Matriz aumentada
    n = nA;    // Tamaño de la matriz
    
    // Eliminación progresiva con pivoteo parcial
    for k = 1:n-1
        kpivot = k; amax = abs(a(k,k));  //pivoteo
        for i = k+1:n
            if abs(a(i, k)) > amax then
                kpivot = i; amax = a(i, k);
            end;
        end;
        temp = a(kpivot, :); a(kpivot, :) = a(k, :); a(k, :) = temp;
        
        for i = k+1:n
            for j = k+1:n+1
                a(i, j) = a(i, j) - a(k, j) * a(i, k) / a(k, k);
            end;
            for j = 1:k        // no hace falta para calcular la solución x
                a(i, j) = 0;  // no hace falta para calcular la solución x
            end              // no hace falta para calcular la solución x
        end;
        it = it + 1;
    end;
    
    // Sustitución regresiva
    x(n) = a(n, n+1) / a(n, n);
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i, k) * x(k);
        end;
        x(i) = (a(i, n+1) - sumk) / a(i, i);
    end;
endfunction

function [p, er] = min_cuad(x, y, n)
    A = mat_min_cuad(x, n);
    [w, a] = gausselimPP((A')*A, (A')*y);
    p = poly(w, "x", "c");
    er = A*w - y;
endfunction
/*
--> x = [0 0.15 0.31 0.5 0.6 0.75]';
--> y = [1 1.001 1.031 1.117 1.223 1.422]';
--> min_cuad(x, y, 3)
 ans  =

                                  2            3
   1.0000684 -0.019602x +0.054392x  +0.9655845x
*/

function A = mat_min_cuad(x, n)
    m = length(x);
    A = ones(m, 1);
    for i = 2:n+1
        A = [A x.^(i-1)];
    end
endfunction
/*
--> x = [0 0.15 0.31 0.5 0.6 0.75]';
--> mat_min_cuad(x, 3)
 ans  =

   1.   0.     0.       0.      
   1.   0.15   0.0225   0.003375
   1.   0.31   0.0961   0.029791
   1.   0.5    0.25     0.125   
   1.   0.6    0.36     0.216   
   1.   0.75   0.5625   0.421875
*/

function z = zeros_chebyshev(n)
    z = zeros(n);
    for i = 1:n
        z(i) = cos((2*i - 1)*%pi / (2*n));
    end
endfunction
/*
--> zeros_chebyshev(4)
 ans  =

   0.9238795
   0.3826834
  -0.3826834
  -0.9238795
*/

function p = cambio_de_variable(g, a, b, n)
    x = zeros_chebyshev(n + 1);
    
    y = zeros(n + 1);
    for i = 1:(n + 1)
        x(i) = ((b + a) + x(i)*(b - a)) / 2;
        y(i) = g(x(i));
    end
    
    p = newton_dd(poly([0], "x"), x, y);
endfunction
/*
--> cambio_de_variable(cos, 0, %pi/2, 3)
 ans  =

                                    2            3
   0.9984416 +0.0319396x -0.6049283x  +0.1142627x
*/

function er = error_poly(p, x, cota)
    n = length(x);
    er = cota/factorial(n);
    for i = 1:n
        er = er * (p - x(i));
    end
endfunction
/*
--> x = [0 0.2 0.4 0.6]';
--> y = [1.0 1.2214 1.4918 1.8221]';
--> error_poly(poly([0], "x"), x, 1.8221)
 ans  =

                         2           3            4
  -0.0036442x +0.0334052x  -0.091105x  +0.0759208x
*/
