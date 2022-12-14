function [x, it] = jacobi(A, b, x0, tol)
    [n, m] = size(A);
    if n <> m then
        error('jacobi - La matriz A debe ser cuadrada');
        abort;
    end;
    it = 0;
    x = x0;
    
    dif = tol + 1;
    
    while dif > tol
        x0 = x;
        
        for i = 1:n
            suma = 0;
            for j = 1:n
                if j <> i
                    suma = suma + A(i, j) * x0(j);
                end
            end
            x(i) = (b(i) - suma) / A(i, i);
        end
        it = it + 1;
        dif = abs(norm(x - x0));
    end
endfunction

function x = jacobi2(A, b, x0, tol)
    [n, m] = size(A);
    x = x0;
    xk = x0;
    
    for i = 1:n
        suma = A(i, 1:n) * xk(1:n);
        suma = suma - A(i, i) * xk(i);
        
        x(i) = (b(i) - suma) / A(i, i);
    end
    
    while abs(norm(x - xk)) > tol
        xk = x;
        
        for i = 1:n
            suma = A(i, 1:n) * xk(1:n);
            suma = suma - A(i, i) * xk(i);
            
            x(i) = (b(i) - suma) / A(i, i)
        end
    end
endfunction

function [x, it] = gauss_seidel(A, b, x0, eps)
    [n, m] = size(A);
    if n <> m then
        error('gauss_seidel - La matriz A debe ser cuadrada');
        abort;
    end;
    
    x = x0;
    it = 0;
    dif = eps + 1;
    
    while dif > eps
        x0 = x;
        
        x(1) = (b(1) - A(1, 2:n) * x0(2:n)) / A(1, 1);
        for i = 2:n-1
            x(i) = (b(i) - A(i, 1:i-1) * x(1:i-1) - A(i, i+1:n) * x0(i+1:n)) / A(i, i);
        end
        x(n) = (b(n) - A(n, 1:n-1) * x(1:n-1)) / A(n, n);
        it = it + 1;
        
        dif = abs(norm(x - x0));
    end
endfunction

function [x, it] = sor(A, b, x0, eps)
    [n, m] = size(A);
    if n <> m then
        error('sor - La matriz A debe ser cuadrada');
        abort;
    end;
    
    T = eye(n, n) - (diag(diag(A))^(-1)) * A;
    p = max(abs(spec(T)));
    w = 2 / (1 + sqrt(1 - p^2));
    
    x = x0;
    it = 0;
    dif = eps + 1;
    
    while dif > eps
        x0 = x;
        
        x(1) = (1 - w) * x0(1) + (w / A(1, 1)) * (b(1) - A(1, 2:n) * x0(2:n));
        for i = 2:n-1
            x(i) = (1 - w) * x0(i) + (w / A(i, i)) * (b(i) - A(i, 1:i-1) * x(1:i-1) - A(i, i+1:n) * x0(i+1:n));
        end
        x(n) = (1 - w) * x0(n) + (w / A(n, n)) * (b(n) - A(n, 1:n-1) * x(1:n-1));
        it = it + 1;
        
        dif = abs(norm(x - x0));
    end
endfunction
