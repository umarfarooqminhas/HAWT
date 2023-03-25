function [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re )

% function used for A_pro, B_pro and C_pro
% interp is a function that allows to find the new values of
% alpha, Cl,Cd, E and Emax for a Reynolds number different from the ones
% given by the data available thanks to a linear interpolation

av=M1(:,1);

n = length(av);
lv = [];
dv = [];
l1 = M1(:,2);
d1 = M1(:,3);
l2 = M2(:,2);
d2 = M2(:,3);

for jkj = 1: n
    lv(jkj) = l1(jkj) + ((l2(jkj)-l1(jkj))/(Re2-Re1))*(Re-Re1);
    dv(jkj) = d1(jkj) + ((d2(jkj)-d1(jkj))/(Re2-Re1))*(Re-Re1);

    %l =
    %lv = [lv;l]
    %j = j+1

end

lv = lv';
dv = dv';

Ev = lv./dv;

Emax = max(Ev);

