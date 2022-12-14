exec('C:\Users\fgozz\OneDrive\Desktop\Facultad\metodos\7\polyChebyshev.sce', -1);
exec('C:\Users\fgozz\OneDrive\Desktop\Facultad\metodos\7\polyNewton.sce', -1);



nodos = rootChebyshev(4);
nodos2 = ((nodos)* (%pi/2) + (%pi/2))/2
w = DD_Newton(nodos2, cos(nodos2))
x = 0 : 0.01: %pi/2
y = cos(x) - horner(w,x)
plot2d(x,y,2)

a=gca()
a.x_location = "origin"
a.y_location = "origin"
 

disp(w)
