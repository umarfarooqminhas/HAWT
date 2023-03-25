function [aopt, lopt, dopt, Eopt] = rad_interp1 (av_1,lv_1,dv_1,Ev_1,av_2,lv_2,dv_2,Ev_2,r ,r1,r2)

% radial interpolation at a generic adimensional radius between two
% profiles

E1max = max(Ev_1);
k1 = find(Ev_1==E1max);
aopt1 = av_1(k1);
lopt1 = lv_1(k1);
dopt1 = dv_1(k1);

E2max = max(Ev_2);
k2 = find(Ev_2==E2max);
aopt2 = av_2(k2);
lopt2 = lv_2(k2);
dopt2 = dv_2(k2);

Eopt = E1max + (E2max-E1max)/(r2-r1)*(r-r1);
aopt = aopt1 + (aopt2-aopt1)/(r2-r1)*(r-r1);
lopt = lopt1 + (lopt2-lopt1)/(r2-r1)*(r-r1);
dopt = dopt1 + (dopt2-dopt1)/(r2-r1)*(r-r1);

