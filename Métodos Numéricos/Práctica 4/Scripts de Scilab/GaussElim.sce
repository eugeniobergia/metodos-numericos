function [x,a] = gausselim(A,b)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana sin pivoteo.  

[nA,mA] = size(A) 
[nb,mb] = size(b)

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
elseif mA<>nb then
    error('gausselim - dimensiones incompatibles entre A y b');
    abort;
end;

a = [A b]; // Matriz aumentada

// Eliminación progresiva
n = nA;
for k=1:n-1
    for i=k+1:n
        a(i,k+1:n+mb) = a(i,k+1:n+mb) - a(k,k+1:n+mb)*a(i,k)/a(k,k);
        a(i,1:k) = 0;
   end;
end;

// Sustitución regresiva
x(n,mb) = 0;
x(n,1:mb) = a(n,n+1:n+mb)/a(n,n);
for i = n-1:-1:1
    sumk = 0;
    sumk = sumk + a(i,i+1:n)*x(i+1:n,1:mb);
    x(i,1:mb) = (a(i,n+1:n+mb)-sumk)/a(i,i);
end;
endfunction

function d = determinante(A)
    [nA,mA] = size(A) 
    if nA<>mA then
        error('gausselim - La matriz A debe ser cuadrada');
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

