function r=raices(p)
    c = coeff(p, 0);
    b = coeff(p, 1);
    a = coeff(p, 2);
    
    d = b**2 - 4*a*c;
    
    if (b > 0) then
        y(1) = (-b - sqrt(d)) / (2*a);
        y(2) = (2*c) / (-b - sqrt(d));
    else if (b < 0) then
            y(1) = (2*c) / (-b + sqrt(d));
            y(2) = (-b + sqrt(d)) / (2*a);
        end
    end
endfunction
