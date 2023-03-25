function [x,k]=newton(x0,nmax,toll,fun,dfun, mol)
% ricerca degli zeri di fun
% xvect = Vettore con tutte le iterate calcolate (xvect(end) = soluzione)

if (nargin == 5)
    mol = 1;     %molteplicità esatta -> si velocizza
end

err = toll + 1;
k = 0;
xv = x0;
xvect = [xv];

while (k< nmax && err> toll)
   fx = fun(xv); 
   dfx = dfun(xv);
   if dfx == 0
      error(' Arresto per azzeramento di dfun');
   else
      xn = xv - mol*(fx/dfx);     
      err = abs(xn-xv);
      xvect = [xvect; xn];
      k = k+1;
      xv = xn;
   end
end

x = xvect(end);