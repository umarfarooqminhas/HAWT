function [av, lv, dv, Ev,Emax] = rad_interp (av_1,lv_1,dv_1,Ev_1,av_2,lv_2,dv_2,Ev_2,r ,r1,r2)

% radial interpolation at a generic adimensional radius between two
% profiles

n = length(av_1);

lv = [];
dv = [];
Ev = [];
av = av_1;

%E1max = max(Ev_1);
%k1 = find(Ev_1==E1max);

%E2max = max(Ev_2);
%k2 = find(Ev_2==E2max);

if length(av_1) == length(av_2)
    
    for jir = 1 : n
        lv(jir) = lv_1(jir)+((lv_2(jir)-lv_1(jir))/(r2-r1))*(r-r1);
        dv(jir) = dv_1(jir)+((dv_2(jir)-dv_1(jir))/(r2-r1))*(r-r1);
        Ev(jir) = Ev_1(jir)+((Ev_2(jir)-Ev_1(jir))/(r2-r1))*(r-r1);

    end

else 
    error('SIR arfa was fucking wrong again')
end

Emax = max(Ev);











