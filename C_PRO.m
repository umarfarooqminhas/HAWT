function [av,lv,dv,Ev,Emax] = C_PRO(Re,C_50K,C_100K,C_200K,C_500K,C_1000K) % C = TIP

% this function provides the new values of alpha, Cl, Cd, E and Emax
% for a given Reynolds number for the hub profile A
% also thanks to the interpolation function interp

av = [];
lv = [];
dv = [];
Ev = [];
Emax=0;

if Re == 50000
    av = C_50K(:,1);
    lv = C_50K(:,2);
    dv = C_50K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 100000
    av = C_100K(:,1);
    lv = C_100K(:,2);
    dv = C_100K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 200000
    av = C_200K(:,1);
    lv = C_200K(:,2);
    dv = C_200K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 500000
    av = C_500K(:,1);
    lv = C_500K(:,2);
    dv = C_500K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 1000000
    av = C_1000K(:,1);
    lv = C_1000K(:,2);
    dv = C_1000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);


elseif Re > 50000 && Re < 100000   % Re btwn 50K and 100K
    M1 = C_50K;
    Re1 = 50000;
    M2 = C_100K;
    Re2 = 100000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re > 100000 && Re < 200000  % Re btwn 100k and 200k
    M1 = C_100K;
    Re1 = 100000;
    M2 = C_200K;
    Re2 = 200000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re > 200000 && Re < 500000  % Re btwn 200k and 500k
    M1 = C_200K;
    Re1 = 200000;
    M2 = C_500K;
    Re2 = 500000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re > 500000 && Re < 1000000  % Re btwn 500k and 1000k
    M1 = C_500K;
    Re1 = 500000;
    M2 = C_1000K;
    Re2 = 1000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re < 50000    % Re lower than the lower reynolds available
    M1 = C_50K;
    Re1 = 50000;
    M2 = C_100K;
    Re2 = 100000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    
elseif Re > 1000000  % Re higher than the higher reynolds available
    av = C_1000K(:,1);
    lv = C_1000K(:,2);
    dv = C_1000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);
end

