function s1 = remonte(A, b)
    n = size(A, 1)
    x(n) = b(n) / A(n, n)
    for i = n - 1:-1:1
        x(i) = (b(i) - A(i, i + 1:n) * x(i + 1:n)) / A(i, i)
    end
    s1 = x
endfunction

function s1 = monte(A, b)
    n = size(A, 1)
    x(n) = 0
    x(1) = b(1) / A(1, 1)
    for i = 2:n
        x(i) = (b(i) - A(i, 1:i - 1) * x(1:i - 1)) / A(i, i)
    end
    s1 = x
endfunction

function [x,a] = gausselim(A,b)
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

function [x,a] = gausselimPP(A,b)
    [nA,mA] = size(A) 
    [nb,mb] = size(b)
    
    if nA<>mA then
        error('gausselimPP - La matriz A debe ser cuadrada');
        abort;
    elseif mA<>nb then
        error('gausselimPP - dimensiones incompatibles entre A y b');
        abort;
    end;
    
    a = [A b]; // Matriz aumentada
    n = nA;    // Tamaño de la matriz
    
    // Eliminación progresiva con pivoteo parcial
    for k=1:n-1
        kpivot = k; amax = abs(a(k,k));  //pivoteo
        for i=k+1:n
            if abs(a(i,k))>amax then
                kpivot = i; amax = a(i,k);
            end;
        end;
        temp = a(kpivot,:); a(kpivot,:) = a(k,:); a(k,:) = temp;
        
        for i=k+1:n
            for j=k+1:n+1
                a(i,j) = a(i,j) - a(k,j)*a(i,k)/a(k,k);
            end;
            for j=1:k        // no hace falta para calcular la solución x
                a(i,j) = 0;  // no hace falta para calcular la solución x
            end              // no hace falta para calcular la solución x
        end;
    end;
    
    // Sustitución regresiva
    x(n) = a(n,n+1)/a(n,n);
    for i = n-1:-1:1
        sumk = 0
        for k=i+1:n
            sumk = sumk + a(i,k)*x(k);
        end;
        x(i) = (a(i,n+1)-sumk)/a(i,i);
    end;
endfunction

function [L, U, P] = lup(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('lup - La matriz A debe ser cuadrada');
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

function [L, U, P] = lup7(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('lup7 - La matriz A debe ser cuadrada');
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
        temp = U(kpivot,k:n); U(kpivot,k:n) = U(k,k:n); U(k,k:n) = temp;
        temp = P(kpivot,:); P(kpivot,:) = P(k,:); P(k,:) = temp;
        temp = L(kpivot,1:k-1); L(kpivot,1:k-1) = L(k,1:k-1); L(k,1:k-1) = temp;
        
        for i=k+1:n
            L(i, k) = U(i,k)/U(k,k);
            U(i,k:n) = U(i,k:n) - U(k,k:n)*L(i, k);
            
            U(i,1:k) = 0;
        end;
    end;
    L = L + eye(n,n);
endfunction

function [L, U] = doolittle(A)
    [nA,mA] = size(A) 
    
    if nA<>mA then
        error('doolittle - La matriz A debe ser cuadrada');
        abort;
    end;
    
    n = nA;
    U = zeros(n, n);
    L = eye(n, n);
    
    U(1,:) = A(1,:);
    L(1:n,1) = A(1:n,1)/U(1,1);
    for i = 2:n
        for j = 2:n
            if i <= j then
                U(i,j) = A(i,j) - L(i,1:i-1)*U(1:i-1,j);
            end
            if i >= j then
                L(i,j) = (A(i,j) - L(i,1:j-1)*U(1:j-1,j))/U(j,j);
            end
        end
    end
endfunction

function U = choleski(A)
    eps = 1.0e-8;
    n = size(A, 1);
    U = zeros(n, n);
    
    t = A(1,1);
    if t <= eps then
        error('choleski - La matriz A no es definida positiva');
        abort;
    end
    U(1,1) = sqrt(t);
    for j = 2:n
        U(1,j) = A(1,j)/U(1,1);
    end
    
    for k = 2:n
        t = A(k,k) - U(1:k-1,k)'*U(1:k-1,k);
        if t <= eps then
            error('choleski - La matriz A no es definida positiva');
            abort;
        end
        U(k,k) = sqrt(t);
        for j = k+1:n
            U(k,j) = (A(k,j) - U(1:k-1,k)'*U(1:k-1,j))/U(k,k);
        end
    end
endfunction

function [Q, R] = my_qr(A)
    [r, c] = size(A);
    if c <> rank(A) then
        error('my_qr - La matriz A no tiene columnas l.i.');
        abort;
    end
    
    Q = zeros(r, c);
    R = zeros(c, c);
    
    v = norm(A(:,1));
    Q(:,1) = A(:,1)/v;
    for i = 1:c
        R(1,i) = A(:,i)' * Q(:,1);
    end
    
    for k = 2:c
        s = zeros(r, 1);
        for i = 1:k-1
            s = s + ((A(:,k)' * Q(:,i)) * Q(:,i));
        end
        v = norm(A(:,k) - s);
        Q(:,k) = (A(:,k) - s)/v;
        for i = k:c
            R(k,i) = A(:,i)' * Q(:,k);
        end
    end
endfunction
