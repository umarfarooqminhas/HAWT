function [av,lv,dv,Ev,Emax] = B_PRO (Re,B_200K,B_500K,B_1000K,B_2000K,B_5000K) 

% this function provides the new values of alpha, Cl, Cd, E and Emax
% for a given Reynolds number for the hub profile A
% also thanks to the interpolation function interp

av = [];
lv = [];
dv = [];
Ev = [];
Emax=0;

if Re == 200000
    av = B_200K(:,1);
    lv = B_200K(:,2);
    dv = B_200K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 500000
    av = B_500K(:,1);
    lv = B_500K(:,2);
    dv = B_500K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 1000000
    av = B_1000K(:,1);
    lv = B_1000K(:,2);
    dv = B_1000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 2000000
    av = B_2000K(:,1);
    lv = B_2000K(:,2);
    dv = B_2000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 5000000
    av = B_5000K(:,1);
    lv = B_5000K(:,2);
    dv = B_5000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re > 200000 && Re < 500000
    M1 = B_200K;
    Re1 = 200000;
    M2 = B_500K;
    Re2 = 500000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re > 500000 && Re < 1000000
    M1 = B_500K;
    Re1 = 500000;
    M2 = B_1000K;
    Re2 = 1000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re > 1000000 && Re < 2000000
    M1 = B_1000K;
    Re1 = 1000000;
    M2 = B_2000K;
    Re2 = 2000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);
    
elseif Re > 2000000 && Re < 5000000
    M1 = B_2000K;
    Re1 = 2000000;
    M2 = B_5000K;
    Re2 = 5000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);
    
    elseif Re < 200000  % Re lower than the lower reynolds available
    M1 = B_200K;
    Re1 = 200000;
    M2 = B_500K;
    Re2 = 500000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );

    elseif Re > 5000000  % Re higher than the higher reynolds available
    av = B_5000K(:,1);
    lv = B_5000K(:,2);
    dv = B_5000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

end











