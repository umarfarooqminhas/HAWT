function [x, niter, err_fin] = fixedpoint (x0,tfun,toll,itmax)

% Authors: R.Sacco, G.Guidoboni and G.Chiaravalli
%
% This function implements the fixed point iteration method
% for finding the zero of a funtion

x = x0;
xold = x;
k = 0;
i = 1;
err = toll + 1;
while ((k < itmax) & (abs(err) > toll))
xnew = tfun(xold);
err = xnew - xold;
x = [x; xnew ];
xold = xnew;
err_fin(i,1)=err; 
k = k+1;
i = i+1;
end

niter = k;

return
