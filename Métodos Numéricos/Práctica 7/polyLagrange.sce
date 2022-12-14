function w = L_i(p,x,i)
    n= length(x);
    w = 1;
    for k= 1: n
        if i <> k then
            w = w* (p - x(k)) /(x(i) -x(k));
        end  
    end
    
endfunction

function w = polyLagrange(p,x,y)
    w = 0;
    n = length(x)
    for i = 1:n
        w = w + y(i)* L_i(p,x,i);
    end
endfunction

x = [0 0.2 0.4 0.6]'
y = [1.0 1.2214 1.4918 1.8221]'

w = polyLagrange(1/3,x,y)
disp(w)
