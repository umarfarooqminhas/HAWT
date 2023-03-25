function [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax)

%function that allows to find the values of alpha, Cl and Cd
% that maximize the efficiency E
% starting for the data extracted by A_pro, B_pro and C_pro

k=1;
%n = length(av);

%if a<av(1)
%    error('alpha is lower than the bottom limit');
%elseif a>av(n)
   % error('alpha is higher than the upper limit');
%elseif a>=av(1) && a<=av(n)


while Emax~= Ev(k) 
        k=k+1;
end

%k =find(Emax==Ev(k)

aopt = av(k);
lopt = lv(k);
dopt = dv(k);

