exec('C:\Users\fgozz\OneDrive\Desktop\Facultad\metodos\4\GaussElimPP.sce', -1); 

function [w,e] = minCuadrados(x,y,p)
    m = length(x)
    A = ones(m,1)
    for j = 2 :p+1
        A = [A, (x'). **(j-1)]
    end
    [w,a] = gaussElimPP(A'*A, A' * y)
    w = poly(w,"x", "coeff") //polinomio
    e = abs(norm(A*w - y,2))
endfunction
