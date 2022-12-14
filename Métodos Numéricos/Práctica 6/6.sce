function [inf, sup] = cotas(A)
    [n, m] = size(A);
    
    centros = diag(A);
    radios = sum(abs(A), 'c') - abs(centros);
    
    inf = min(centros - radios);
    sup = max(centros + radios);
endfunction

function c_p = c_poly(A)
    [n, m] = size(A);
    if n <> m then
        error('c_poly - La matriz A debe ser cuadrada');
        abort;
    end;
    p = poly([0 1], 'x', 'c');
    c_p = det(A - eye(n, n) * p);
endfunction

function circ(r, x, y)
    xarc(x - r, y + r, 2 * r, 2 * r, 0, 360 * 64);
endfunction

function gres(A)
    [n, m] = size(A);
    
    centros = diag(A);
    radios = sum(abs(A), 'c') - abs(centros);
    
    mx = round(min(centros - radios) - 1);
    Mx = round(max(centros + radios) + 1);
    
    my = round(min(-radios) - 1);
    My = -my;
    
    rectangulo = [mx, my, Mx, My];
    
    s = spec(A);
    plot2d(real(s), imag(s), -1, "031", "", rectangulo);
    // replot(rectangulo);
    xgrid();
    for i = 1:n
        circ(radios(i), centros(i), 0);
    end
endfunction

function e = potencia(A, z0, max_it)
    [n, m] = size(A);
    if n <> m then
        error('potencia - La matriz A debe ser cuadrada');
        abort;
    end;
    
    z = z0;
    w = z0;
    it = 0;
    
    while it < max_it
        z0 = z;
        w = A * z0;
        
        z = w / norm(w, %inf);
        
        it = it + 1;
    end
    
    k = 0;
    flag = 0;
    for i = 1:n
        if w(i) <> 0 then
            k = i;
            flag = 1;
            break;
        end;
    end;
    if flag == 0 then
        error('potencia - No es posible aproximar el radio espectral');
        abort;
    end
    
    e = w(k) / z0(k);
endfunction

function [e, z, it] = potencia2(A, z0, eps, max_it)
    [n, m] = size(A);
    if n <> m then
        error('potencia2 - La matriz A debe ser cuadrada');
        abort;
    end;
    
    z = z0;
    w = z0;
    dif = eps + 1;
    it = 0;
    
    while dif > eps && it < max_it
        z0 = z;
        w = A*z0;
        
        z = w / norm(w, %inf);
        
        [ma, k] = max(abs(w));
        e = w(k) / z0(k);
        
        dif = abs(norm(z - z0));
        it = it + 1;
    end
endfunction
