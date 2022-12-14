function s1 = remonte(A, b)
    n = size(A, 1);
    x(n) = b(n) / A(n, n);
    for i = n - 1:-1:1
        x(i) = (b(i) - A(i, i + 1:n) * x(i + 1:n)) / A(i, i);
    end
    s1 = x;
endfunction

function s1 = monte(A, b)
    n = size(A, 1);
    x(n) = 0;
    x(1) = b(1) / A(1, 1);
    for i = 2:n
        x(i) = (b(i) - A(i, 1:i - 1) * x(1:i - 1)) / A(i, i);
    end
    s1 = x;
endfunction

function [x, a, it] = gausselim(A, b)
    [nA, mA] = size(A);
    [nb, mb] = size(b);
    
    if nA <> mA then
        error('gausselim - La matriz A debe ser cuadrada');
        abort;
    elseif mA <> nb then
        error('gausselim - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    it = 0;
    a = [A b]; // Matriz aumentada
    
    // Eliminación progresiva
    n = nA;
    for k=1:n-1
        for i=k+1:n
            a(i, k+1:n+mb) = a(i, k+1:n+mb) - a(k, k+1:n+mb) * a(i, k) / a(k, k);
            a(i, 1:k) = 0;
            it = it + 1;
       end;
    end;
    
    // Sustitución regresiva
    x(n, mb) = 0;
    x(n, 1:mb) = a(n, n+1:n+mb) / a(n, n);
    for i = n-1:-1:1
        sumk = 0;
        sumk = sumk + a(i, i+1:n) * x(i+1:n, 1:mb);
        x(i, 1:mb) = (a(i, n+1:n+mb) - sumk) / a(i, i);
    end;
endfunction
/*
Ejemplo:

A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1];
b = [4 1 -3 4]';
[x, a, it] = gausselim(A, b)
 it  = 

   6.

 a  = 

   1.   1.   0.   3.    4. 
   0.  -1.  -1.  -5.   -7. 
   0.   0.   3.   13.   13.
   0.   0.   0.  -13.  -13.

 x  = 

  -1.
   2.
   0.
   1.
*/

function d = determinante(A)
    [nA, mA] = size(A) 
    if nA <> mA then
        error('determinante - La matriz A debe ser cuadrada');
        abort;
    end;
    
    n = nA;
    a = A;
    for k = 1:n-1
        for i = k+1:n
            a(i, k+1:n) = a(i, k+1:n) - a(k, k+1:n) * a(i, k) / a(k, k);
            a(i, 1:k) = 0;
       end;
    end;
    
    d = 1;
    for i = 1:n
        d = d * a(i, i);
    end
endfunction
/*
Ejemplo:

A = [1 1 0 3; 2 1 -1 1; 3 -1 -1 2; -1 2 3 -1];
d = determinante(A)
 d  = 

   39.
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
/*
Ejemplo:

A = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b = [-8 -20 -2 4]';
[x, a, it] = gausselimPP(A, b)
 it  = 

   3.

 a  = 

   2.  -2.   3.   -3.   -20.
   0.   2.  -0.5   1.5   8. 
   0.   0.   2.5   4.5   14.
   0.   0.   0.   -0.4  -0.8

 x  = 

  -7.
   3.
   2.
   2.
*/

function [L, U, P] = lup(A)
    [nA, mA] = size(A) 
    
    if nA <> mA then
        error('lup - La matriz A debe ser cuadrada');
        abort;
    end;
    
    U = A; // Matriz aumentada
    n = nA;    // Tamaño de la matriz
    P = eye(n, n);
    L = zeros(n, n);
    
    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(U(k, k));  //pivoteo
        for i=k+1:n
            if abs(U(i, k)) > amax then
                kpivot = i; amax = U(i, k);
            end;
        end;
        temp = U(kpivot,:); U(kpivot,:) = U(k,:); U(k,:) = temp;
        temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp;
        temp = L(kpivot,:); L(kpivot,:) = L(k,:); L(k,:) = temp;
        
        for i=k+1:n
            m = U(i, k) / U(k, k);
            U(i, k+1:n) = U(i, k+1:n) - U(k, k+1:n)*m;
            L(i, k) = m;
            
            U(i, 1:k) = 0;
        end;
    end;
    L = L + eye(n, n);
endfunction
/*
Ejemplo:

A = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b = [-8 -20 -2 4]';
[L, U, P] = lup(A)
 P  = 

   0.   1.   0.   0.
   0.   0.   1.   0.
   0.   0.   0.   1.
   1.   0.   0.   0.

 U  = 

   2.  -2.   3.   -3. 
   0.   2.  -0.5   1.5
   0.   0.   2.5   4.5
   0.   0.   0.   -0.4

 L  = 

   1.    0.   0.    0.
   0.5   1.   0.    0.
   0.5   0.   1.    0.
   0.5   0.   0.2   1.
*/

function [L, U, P] = lup7(A)
    [nA, mA] = size(A) 
    
    if nA <> mA then
        error('lup7 - La matriz A debe ser cuadrada');
        abort;
    end;
    
    U = A; // Matriz aumentada
    n = nA;    // Tamaño de la matriz
    P = eye(n, n);
    L = zeros(n, n);
    
    // Eliminación progresiva con pivoteo parcial
    for k = 1:n-1
        kpivot = k; amax = abs(U(k, k));  //pivoteo
        for i = k+1:n
            if abs(U(i, k))>amax then
                kpivot = i; amax = U(i, k);
            end;
        end;
        temp = U(kpivot, k:n); U(kpivot, k:n) = U(k, k:n); U(k, k:n) = temp;
        temp = P(kpivot, :); P(kpivot, :) = P(k, :); P(k, :) = temp;
        temp = L(kpivot, 1:k-1); L(kpivot, 1:k-1) = L(k, 1:k-1); L(k, 1:k-1) = temp;
        
        for i= k+1:n
            L(i, k) = U(i, k) / U(k, k);
            U(i, k:n) = U(i, k:n) - U(k, k:n) * L(i, k);
            
            U(i, 1:k) = 0;
        end;
    end;
    L = L + eye(n, n);
endfunction
/*
Ejemplo:

A = [1 -1 2 -1; 2 -2 3 -3; 1 1 1 0; 1 -1 4 3]
b = [-8 -20 -2 4]';
[L, U, P] = lup7(A)
 P  = 

   0.   1.   0.   0.
   0.   0.   1.   0.
   0.   0.   0.   1.
   1.   0.   0.   0.

 U  = 

   2.  -2.   3.   -3. 
   0.   2.  -0.5   1.5
   0.   0.   2.5   4.5
   0.   0.   0.   -0.4

 L  = 

   1.    0.   0.    0.
   0.5   1.   0.    0.
   0.5   0.   1.    0.
   0.5   0.   0.2   1.
*/

function [L, U] = doolittle(A)
    [nA, mA] = size(A) 
    
    if nA <> mA then
        error('doolittle - La matriz A debe ser cuadrada');
        abort;
    end;
    
    n = nA;
    U = zeros(n, n);
    L = eye(n, n);
    
    U(1, :) = A(1, :);
    L(1: n,1) = A(1:n, 1) / U(1, 1);
    for i = 2:n
        for j = 2:n
            if i <= j then
                U(i, j) = A(i, j) - L(i, 1:i-1) * U(1:i-1, j);
            end
            if i >= j then
                L(i, j) = (A(i, j) - L(i, 1:j-1) * U(1:j-1, j)) / U(j, j);
            end
        end
    end
endfunction
/*
Ejemplo:

A = [1 2 3 4; 1 4 9 16; 1 8 27 64; 1 16 81 256];
[L, U] = doolittle(A)
 U  = 

   1.   2.   3.   4. 
   0.   2.   6.   12.
   0.   0.   6.   24.
   0.   0.   0.   24.

 L  = 

   1.   0.   0.   0.
   1.   1.   0.   0.
   1.   3.   1.   0.
   1.   7.   6.   1.
*/

function U = choleski(A)
    eps = 1.0e-8;
    n = size(A, 1);
    U = zeros(n, n);
    
    t = A(1, 1);
    if t <= eps then
        error('choleski - La matriz A no es definida positiva');
        abort;
    end
    U(1, 1) = sqrt(t);
    for j = 2:n
        U(1, j) = A(1, j) / U(1, 1);
    end
    
    for k = 2:n
        t = A(k, k) - U(1:k-1, k)' * U(1:k-1, k);
        if t <= eps then
            error('choleski - La matriz A no es definida positiva');
            abort;
        end
        U(k, k) = sqrt(t);
        for j = k+1:n
            U(k, j) = (A(k, j) - U(1:k-1, k)' * U(1:k-1, j)) / U(k, k);
        end
    end
endfunction
/*
Ejemplo:

A = [16 -12 8 -16; -12 18 -6 9; 8 -6 5 -10; -16 9 -10 46];
U = choleski(A)
 U  = 

   4.  -3.   2.  -4.
   0.   3.   0.  -1.
   0.   0.   1.  -2.
   0.   0.   0.   5.
*/

function [Q, R] = my_qr(A)
    [r, c] = size(A);
    if c <> rank(A) then
        error('my_qr - La matriz A no tiene columnas l.i.');
        abort;
    end
    
    for k = 1:c
        s = zeros(r, 1);
        for i = 1:k-1
            s = s + ((A(:, k)' * Q(:, i)) * Q(:, i));
        end
        v = norm(A(:, k) - s);
        Q(:, k) = (A(:, k) - s) / v;
        for i = k:c
            R(k, i) = A(:, i)' * Q(:, k);
        end
    end
endfunction
/*
Ejemplo:

A = [16 -12 8; -12 18 -6; 8 -6 8];
[Q, R] = my_qr(A)
 R  = 

   21.540659  -21.169269   12.255892
   0.          7.4740932   0.9965458
   0.          0.          3.5777088

 Q  = 

   0.7427814   0.4982729  -0.4472136
  -0.557086    0.8304548   2.483D-16
   0.3713907   0.2491364   0.8944272
*/