function integral = gauss_quad(a, b, fun, Mh)

% integral = gauss_quad(a, b, fun, Mh)
% Note: fun is a handle function whose syntax is fun = @(x) f(x)

h = (b-a)/Mh;

ab_v= [a:h:b];

integral = 0;

for K = 1:Mh

    xa = a + (K-1)*h;
	xb = xa + h;
	xc = (xa + xb)/2;

	xg1 = xc - h/(2*sqrt(3));
	xg2 = xc + h/(2*sqrt(3));
	
    int_K = h/2*(fun(find_tol(ab_v,xg1,0.1)) + fun(find_tol(ab_v,xg1,0.1)));
    integral = integral + int_K;
end
%
return
