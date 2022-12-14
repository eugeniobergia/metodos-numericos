function [L, U, P] = lup(A)
// Esta función obtiene la solución del sistema de ecuaciones lineales A*x=b, 
// dada la matriz de coeficientes A y el vector b.
// La función implementa el método de Eliminación Gaussiana con pivoteo parcial.

[nA,mA] = size(A) 

if nA<>mA then
    error('gausselim - La matriz A debe ser cuadrada');
    abort;
end;

U = A; // Matriz aumentada
n = nA;    // Tamaño de la matriz
P = eye(n, n);
L = zeros(n, n);

// Eliminación progresiva con pivoteo parcial
for k=1:n-1
    kpivot = k; amax = abs(U(k,k));  //pivoteo
    for i=k+1:n
        if abs(U(i,k))>amax then
            kpivot = i; amax = U(i,k);
        end;
    end;
    temp = U(kpivot,:); U(kpivot,:) = U(k,:); U(k,:) = temp;
    temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp;
    temp = L(kpivot,:); L(kpivot,:) = L(k,:); L(k,:) = temp;
    
    for i=k+1:n
        m = U(i,k)/U(k,k);
        U(i,k+1:n) = U(i,k+1:n) - U(k,k+1:n)*m;
        L(i, k) = m;
        
        U(i,1:k) = 0;
    end;
end;
L = L + eye(n,n);

endfunction
