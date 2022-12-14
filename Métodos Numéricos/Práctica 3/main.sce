funcprot(0);

function ret = ej9(x0)
    x = x0(1);
    y = x0(2);
    f1 = 1 + x**2 -y**2 + %e**x * cos(y);
    f2 = 2 * x * y + %e**x * sin(y);
    
    ret = [f1; f2]
endfunction

function ret = ej11(x0)
    x = x0(1);
    y = x0(2);
    f1 = 2 + 4 * x * %e**(2 * x**2 + y**2);
    f2 = 6 * y + 2 * y * %e**(2 * x**2 + y**2);
    
    ret = [f1; f2]
endfunction

function ret = biseccion(f, a, b, tol)
    c = (a + b)/2;
    while b - c > tol
        if f(c)*f(b) <= 0 then
            a = c;
        else
            b = c;
        end
        c = (a + b)/2;
    end
    ret = c;
endfunction

function ret = biseccion2(f, a, b, it, tol)
    c = (a + b)/2;
    i = 1;
    while i < it && b - c > tol
        if f(c)*f(b) <= 0 then
            a = c;
        else
            b = c;
        end
        c = (a + b)/2;
        i = i + 1;
    end
    ret = c;
endfunction

function ret = newton(f, x0, it, tol)
    xi = x0 - (numderivative(f, x0)**(-1)) * f(x0);
    i = 0;
    while i < it && abs(xi - x0) > tol
        x0 = xi;
        xi = x0 - (numderivative(f, x0)**(-1)) * f(x0);
        i = i + 1;
    end
    ret = xi;
endfunction

function ret = secante(f, x0, x1, tol)
    xi = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)));
    while abs(xi - x1) > tol
        x0 = x1;
        x1 = xi;
        xi = x1 - f(x1) * ((x1 - x0) / (f(x1) - f(x0)));
    end
    ret = xi;
endfunction

function ret = regula_falsi(f, a, b, tol)
    c = b - f(b) * ((b - a) / (f(b) - f(a)));
    while b - c > tol
        if f(c)*f(b) <= 0 then
            a = c;
        else
            b = c;
        end
        c = b - f(b) * ((b - a) / (f(b) - f(a)));
    end
    ret = c;
endfunction

function ret = punto_fijo(f, x0, tol)
    xi = f(x0);
    while abs(xi - x0) > tol
        x0 = xi;
        xi = f(x0);
    end
    ret = xi;
endfunction
