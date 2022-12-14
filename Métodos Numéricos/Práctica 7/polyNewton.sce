function w=DD(x,y)
    n = length(x)
    if n==2 then
        w = ( y(2)-y(1) ) / (x(2)- x(1))
    else
        w = ( DD(x(2:n), y(2:n)) - DD( x(1:n-1), y(1:n-1)) ) / (x(n) - x(1))
    end
    
endfunction


function w = polyNewton(p,x,y)
    w = 0;
    n = length(x);
    for j = n: -1 :2
        w = (w+ DD(x(1:j), y(1:j))) * (p-x(j-1)) ;
    end
    w = w + y(1);
    
endfunction


function w = err(p,x,cot)
    n = length(x)
    w = cot / factorial(n)
    for i = 1:n
        w = w*(p - x(i))
    end
endfunction

function w = DD_Newton(x,y)
    // Entrada: x,y = vectores puntos de interpolación (x,y)
    // Salida: w = polinomio de diferencias divididas de Newton
    w = 0
    s = poly(0,"x")
    n = length(x)
    for j=n:-1:2
        w = w*(s-x(j-1)) + DD(x(1:j),y(1:j))*(s-x(j-1))
    end
    w = w + y(1)
endfunction


//x = [0 0.2 0.4 0.6]'
//y = [1.0 1.2214 1.4918 1.8221]'
//
//w = polyNewton(1/3,x,y)
//disp(w)
