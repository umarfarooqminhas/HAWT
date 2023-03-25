function [av,lv,dv,Ev,Emax] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K) % A = ROOT

% this function provides the new values of alpha, Cl, Cd, E and Emax
% for a given Reynolds number for the hub profile A
% also thanks to the interpolation function interp

av = [];
lv = [];
dv = [];
Ev = [];
Emax=0;

if Re == 200000
    av = A_200K(:,1);
    lv = A_200K(:,2);
    dv = A_200K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 500000
    av = A_500K(:,1);
    lv = A_500K(:,2);
    dv = A_500K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 1000000
    av = A_1000K(:,1);
    lv = A_1000K(:,2);
    dv = A_1000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 2000000
    av = A_2000K(:,1);
    lv = A_2000K(:,2);
    dv = A_2000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);

elseif Re == 5000000
    av = A_5000K(:,1);
    lv = A_5000K(:,2);
    dv = A_5000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);


elseif Re > 200000 && Re < 500000   % Re btwn 20K and 50K
    M1 = A_200K;
    Re1 = 200000;
    M2 = A_500K;
    Re2 = 500000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    
elseif Re > 500000 && Re < 1000000  % Re btwn 50k and 100k
    M1 = A_500K;
    Re1 = 500000;
    M2 = A_1000K;
    Re2 = 1000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
   %[a, l, d] = finder(a,av,lv,dv);
   %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re > 1000000 && Re < 2000000  % Re btwn 100k and 200k
    M1 = A_1000K;
    Re1 = 1000000;
    M2 = A_2000K;
    Re2 = 2000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[a, l, d] = finder(a,av,lv,dv);
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re > 2000000 && Re < 5000000  % Re btwn 200k and 500k
    M1 = A_2000K;
    Re1 = 2000000;
    M2 = A_5000K;
    Re2 = 5000000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );
    %[a, l, d] = finder(a,av,lv,dv);
    %[aopt, lopt, dopt] = finder(av, lv, dv, Emax, E);

elseif Re < 200000    % Re lower than the lower reynolds available
    M1 = A_200K;
    Re1 = 200000;
    M2 = A_500K;
    Re2 = 500000;
    [av, lv, dv, Ev, Emax] = interp (M1 ,M2 ,Re1, Re2 ,Re );

elseif Re > 5000000    % Re higher than the higher reynolds available
    av = A_5000K(:,1);
    lv = A_5000K(:,2);
    dv = A_5000K(:,3);
    Ev = lv./dv;
    Emax = max(Ev);
end









