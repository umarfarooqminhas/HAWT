%% OFF DESIGN - 3 BLADED ROTOR HAWT

V_A = 3;
V_cutin = V_A;

V_E = 27;
V_cutout = V_E;

u_t_max = 80; 
R = 55;
omega_max_rps = u_t_max/R; %[rps]
omega_max = omega_max_rps *60; %rpm
V0_on=7.3;
V_D = V0_on * 1.7;
V_rated = V_D;

lambda_max = 6.9; % lam = 6.9

V_C = u_t_max/lambda_max;
V_cpmax = V_C;

V0_v = linspace(V_A,V_E, 50);
V0_v = V0_v';

cp_v=[];
P_voff=[];

bc_on = bc_v; %[°]

j_A = find_tol(V0_v,V_A,0.25);
j_C = find_tol(V0_v,V_C,0.25);
j_D = find_tol(V0_v,V_D,0.25);
j_E = length(V0_v);

%% VA < V0 < VC

lam = lambda_max;

for j_ac_off = j_A : j_C

    V0 = V0_v(j_ac_off);
    omega = lam*V0/R; %rps
   
    [cp, P_AC] = off_AC (bc_on, lam, V0, a_corrected, ap_corrected, cv_final);
    cp_v(j_ac_off)= cp;
    P_voff(j_ac_off) = P_AC;
    omega_v(j_ac_off) = omega;

end

%% VC < V0 < VD

for j_cd_off = j_C+1: j_D
    
    V0 = V0_v(j_cd_off);
    omega = omega_max_rps; 
    lam = omega*R/V0;
    D_bc = 0;
    [cp, P_CD] = off_CD (bc_on, D_bc , lam, V0, a_corrected, ap_corrected, cv_final);
    
    cp_v(j_cd_off)= cp;
    P_voff(j_cd_off) = P_CD;

end

%% VD < V0 <VE

omega = omega_max_rps; 
P_rated = P_voff(j_D);

for j_de_off = j_D+1 : j_E
    
    V0 = V0_v(j_de_off);
    V0
    lam = omega*R/V0;
    D_bc = 0;
    [cp, P_DE] = off_DE (bc_on, D_bc , lam, V0, a_corrected, ap_corrected, cv_final, P_rated);
    cp_v(j_de_off)= cp;
    P_voff(j_de_off) = P_DE;

end
   

%% Plots
figure;
plot(V0_v(j_A:j_E)',cp_v(j_A:j_E),'LineWidth',2)
xlabel('V [m/s]');
ylabel('Cp[-]')
title('Power Coefficient Vs Wind Velocity')
%plot(V0_v',cp_v)
figure;
plot(V0_v(j_A:j_E),P_voff(j_A:j_E),'LineWidth',2);
xlabel('V [m/s]');
ylabel('P[W]')
title('Power Vs Wind Velocity')








