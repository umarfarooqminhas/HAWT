%% ON DESIGN - 3 BLADED ROTOR HAWT
clear 
close all
clc

%% DATA EXTRACTION
s=0.05;
a_first = -5;
a_end = 15;

%ROOT

load A_sp4721_200K.txt;
load A_sp4721_500K.txt;
load A_sp4721_1000K.txt;
load A_sp4721_2000K.txt;
load A_sp4721_5000K.txt;

[A_200K]  = palas(A_sp4721_200K,s, a_first, a_end);  %A_sp4721_200K
[A_500K]  = palas(A_sp4721_500K,s, a_first, a_end);  %A_sp4721_500K
[A_1000K] = palas(A_sp4721_1000K,s, a_first, a_end); %A_sp4721_1000K
[A_2000K] = palas(A_sp4721_2000K,s, a_first, a_end); %A_sp4721_2000K
[A_5000K] = palas(A_sp4721_5000K,s, a_first, a_end); %A_sp4721_5000K

%PRIMARY

load B_sp4621_200K.txt;
load B_sp4621_500K.txt;
load B_sp4621_1000K.txt;
load B_sp4621_2000K.txt;
load B_sp4621_5000K.txt;

[B_200K]  = palas(B_sp4621_200K,s, a_first, a_end);  %B_sp4621_200K
[B_500K]  = palas(B_sp4621_500K,s, a_first, a_end);  %B_sp4621_500K
[B_1000K] = palas(B_sp4621_1000K,s, a_first, a_end); %B_sp4621_1000K
[B_2000K] = palas(B_sp4621_2000K,s, a_first, a_end); %B_sp4621_2000K
[B_5000K] = palas(B_sp4621_5000K,s, a_first, a_end); %B_sp4621_5000K

%TIP

load C_S828_50K.txt;
load C_S828_100K.txt;
load C_S828_200K.txt;
load C_S828_500K.txt;
load C_S828_1000K.txt;

[C_50K]   = palas(C_S828_50K,s, a_first, a_end);   %C_S828_50K
[C_100K]  = palas(C_S828_100K,s, a_first, a_end);  %C_S828_100K
[C_200K]  = palas(C_S828_200K,s, a_first, a_end);  %C_S828_200K
[C_500K]  = palas(C_S828_500K,s, a_first, a_end);  %C_S828_500K
[C_1000K] = palas(C_S828_1000K,s, a_first, a_end); %C_S828_1000K

%% DATA & assumptions

D=110;                    % DIAMETER	                [m]
R = D/2;                  % Radius                      [m]
Nb=3;	                  % Number of Blades            [-]		
mi=1.81*10^(-5);          % Dynamic viscosity           [Kg/ms]
ro=1.22;                  % Rho = density               [kg/m3]
v=mi / ro;                % Kinematic Viscosity         [m2/s]
z = 100;                  % altidude	                [m]
V0 = 7.30;                % Wind velocity (undisturbed) [m/s]
lam = 6.90;               % lambda = tip speed ratio    [-]
rot_m =	V0*lam/R*60;      % Rotational speed            [rpm]
rot_s = rot_m/60;         % Rotational speed            [rps]
rot_rad = rot_s*2*pi;     % Rotational speed            [rad/s]
rot_gra = rot_rad*180/pi; % Rotational speed            [°/s]

omega=rot_m; %it may be changed according to the needs


%% Glauert Distribution            %0.0001  % .000001 .000000001  -0

a_gv = [0.260: 0.000000001 : 0.3328]'; % max = 0.3335 otherwise complex numbers
ap_gv = (1-3*a_gv)./(4*a_gv-1);
xv = (a_gv.*(1-a_gv)./(ap_gv.*(1.+ap_gv))).^0.5;
fi_rad = atan(ap_gv.*xv./a_gv);
fi_deg = fi_rad.*180./pi;
w_infv = V0.*(1-a_gv)./sin(fi_rad);

r_ad = xv./lam; % r/R = r_adimensional [-]
fi_v=fi_rad;

rh = 0.10; % adimensional radius ratio - hub
ra = 0.20; % ROOT PROFILE = A
rr = 0.40; % adimensional radius ratio - root
rb = 0.50; % PRIMARY PROFILE = B
rp = 0.80; % adimensional radius ratio - primary
rc = 0.86; % TIP PROFILE = C
rt = 0.95; % adimensional radius ratio - tip
rf = r_ad(end); % adimensional radius ratio - final

R_h = rh * R ; % R_hub
R_r = rr * R ; % R_root
R_p = rp * R ; % R_primary
R_t = rt * R ; % R_tip
R_f = rf * R ; % R_final

%% tolerances selection

drad = 0.01;          % 0.01 0.001

tol = 0.01;           % 0.01

tol_h = 0.01;         % 0.01
tol_a = 0.01;          % 0.1
tol_b = 0.01;    %0.2 %0.3 no => %0.25 
tol_c = 0.01;    %0.3  %0.003 %000023   0000025 000003

%% optimal chord distribution 

%% EFFECT OF PROFILE A

r_ad_A = [rh:drad:ra]';
c=5; %first guess of c => c @ hub [m]

err= tol +1;

while err > tol
    x = rh*lam;
    [kh] = find_tol (xv,x,tol_h);
    a = a_gv(kh);
    ap= ap_gv(kh);
    winf = w_infv(kh);
    fi = fi_rad(kh);
    Re=c*winf/v;
    
    [av,lv,dv,Ev,Emax] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
    
    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);

    F = (2/pi)*acos(exp(Nb*(rh-1)/(2*rh*sin(fi))));

    c_n=(4*F/lopt)*(2*pi*rh*R/Nb)*(a/(1-a))*sin(fi)*tan(fi);
    
    err = abs(c - c_n);

    c=c_n;
end

cv(1) = c;
cl_vect(1) = lopt;
a_v(1)=a;

for j = 2:length(r_ad_A)
    
    co= cv(j-1);
    err=err+1;
    
    while err > tol
    
    x = r_ad_A(j)*lam;
    [ka] = find_tol (xv,x,tol_a);
    a = a_gv(ka);
    %a
    ap= ap_gv(ka);
    winf = w_infv(ka);
    fi = fi_rad(ka);
    Re=co*winf/v;
    
    [av_A,lv_A,dv_A,Ev_A,Emax_A] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
  
    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
   
    F = (2/pi)*acos(exp(Nb*(r_ad_A(j)-1)/(2*r_ad_A(j)*sin(fi))));
    c_n=(4*F/lopt)*(2*pi*r_ad_A(j)*R/Nb)*(a/(1-a))*sin(fi)*tan(fi);
    
    err = abs(co - c_n);
    
    co=c_n;

    end

    cv(j) = co;
    cl_vect(j) = lopt;
    avv(j)=a;
    a_v(j)=a;
   
end

ja=length(cv);

%% EFFECT OF A & B

r_ad_AB = [ra+drad:drad:rb]';

for j = 1:length(r_ad_AB)
    
    co= cv(j);
    err=err+1;
    
    while err > tol
    
    x = r_ad_AB(j)*lam;
    [kab] = find_tol (xv,x,tol_a);
    a = a_gv(kab);
    ap= ap_gv(kab);
    winf = w_infv(kab);
    fi = fi_rad(kab);
    Re=co*winf/v;
     
    [av_A,lv_A,dv_A,Ev_A,Emax_A] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
    [av_B,lv_B,dv_B,Ev_B,Emax_B] = B_PRO(Re,B_200K,B_500K,B_1000K,B_2000K,B_5000K);
    
    [av, lv, dv, Ev,Emax] = rad_interp(av_A,lv_A,dv_A,Ev_A,av_B,lv_B,dv_B,Ev_B,r_ad_AB(j),ra,rb);
   
    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
    
    F = (2/pi)*acos(exp(Nb*(r_ad_AB(j)-1)/(2*r_ad_AB(j)*sin(fi))));
    c_n=(4*F/lopt)*(2*pi*r_ad_AB(j)*R/Nb)*(a/(1-a))*sin(fi)*tan(fi);
    
    err = abs(co - c_n);
    co=c_n;

    end

    cv(j+ja) = co;
    cl_vect(j+ja) = lopt;
    avv(j+ja)=a;
    a_v(j+ja)=a;
    
end

jab=length(cv);

%% optimal chord distribution - PRIMARY

r_ad_BC = [rb+drad:drad:rc]';

for j = 1:length(r_ad_BC)
    
    co= cv(j-1+jab);
    err=err+1;
    
    while err > tol
    
    x = r_ad_BC(j)*lam;
    [kbc] = find_tol (xv,x,tol_b); 
    a = a_gv(kbc);
    ap= ap_gv(kbc);
    winf = w_infv(kbc);
    fi = fi_rad(kbc);
    Re=co*winf/v;
    
    [av_B,lv_B,dv_B,Ev_B,Emax_B] = B_PRO(Re,B_200K,B_500K,B_1000K,B_2000K,B_5000K);
    [av_C,lv_C,dv_C,Ev_C,Emax_C] = C_PRO(Re,C_50K,C_100K,C_200K,C_500K,C_1000K);
    
    [av, lv, dv, Ev,Emax] = rad_interp(av_B,lv_B,dv_B,Ev_B,av_C,lv_C,dv_C,Ev_C,r_ad_BC(j),rb,rc);
   
    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
    
    F = (2/pi)*acos(exp(Nb*(r_ad_BC(j)-1)/(2*r_ad_BC(j)*sin(fi))));
    
    c_n=(4*F/lopt)*(2*pi*r_ad_BC(j)*R/Nb)*(a/(1-a))*sin(fi)*tan(fi);

    err = abs(co - c_n);
    
    co=c_n;
    
    end

    cv(jab+j) = co;
    cl_vect(jab+j) = lopt;
    a_v(j+jab)=a;

end

jbc=length(cv);

%% optimal chord distribution - TIP

r_ad_C = [rc+drad:drad:rt]';

for j = 1:length(r_ad_C)

    co= cv(j-1+jbc);
    err=err+1;
    
    while err > tol
    
    x = r_ad_C(j)*lam;
    [kc] = find_tol (xv,x,tol_c); %0.3   %0.003 %000023
    a = a_gv(kc);
    ap= ap_gv(kc);
    winf = w_infv(kc);
    fi = fi_rad(kc);
    Re=co*winf/v;
    
    [av_C,lv_C,dv_C,Ev_C,Emax_C] = C_PRO(Re,C_50K,C_100K,C_200K,C_500K,C_1000K);
    
    [aopt, lopt, dopt] = finder(av_C, lv_C, dv_C, Ev_C, Emax_C);
    F = (2/pi)*acos(exp(Nb*(r_ad_C(j)-1)/(2*r_ad_C(j)*sin(fi))));
    c_n=(4*F/lopt)*(2*pi*r_ad_C(j)*R/Nb)*(a/(1-a))*sin(fi)*tan(fi);

    err = abs(co - c_n);
    
    co=c_n;
    
    end

    cv(jbc+j) = co;
    cl_vect(jbc+j) = lopt;
    a_v(j+jbc)=a;

end
cv = cv';

%%
figure;
ncv = length(cv);
r_ad_tot=[rh:drad:rt]';
r_ad_interp = interpft(r_ad_tot, ncv);
plot(r_ad_interp,cv,'r','LineWidth',2);
title('Optimal Chord Distribution');
xlabel('c/R');
ylabel('r/R');

%% linearization
cvrad= cv./R;
%figure
%plot(r_ad_t, cvrad);
title('cv / R & r_a_d_,_t');

%k = find(r_ad_t == 0.4);
%cvrad_new=[cvrad(k):0.001:cvrad]

cvrad_h = 0.12;

kb = find(r_ad_tot == rb);
cvrad_b = cvrad(kb);

kc = find(r_ad_tot == rc);
cvrad_c = cvrad(kc);

kt = find(r_ad_tot == rt);
cvrad_t = cvrad(kt);

r_ad_lin1 = [rh:drad:rb]';
r_ad_lin2 = [rb+drad:drad:rc]';
r_ad_lin3 = [rc+drad:drad:rt]';

cvrad_lin1 = @(r_ad) ((cvrad_b-cvrad_h)/(rb-rh)).*(r_ad-rh) + cvrad_h;
cvrad_lin2 = @(r_ad) ((cvrad_c-cvrad_b)/(rc-rb)).*(r_ad-rb) + cvrad_b;
cvrad_lin3 = @(r_ad) ((cvrad_t-cvrad_c)/(rt-rc)).*(r_ad-rc) + cvrad_c;

cv_ad_lin1=cvrad_lin1(r_ad_lin1);
cv_ad_lin2=cvrad_lin2(r_ad_lin2);
cv_ad_lin3=cvrad_lin3(r_ad_lin3);

cv_ad_lin = [cv_ad_lin1;cv_ad_lin2;cv_ad_lin3];
r_ad_lin  = [r_ad_lin1;r_ad_lin2;r_ad_lin3 ];

figure
plot(r_ad_interp,cv./R,'k','LineWidth',2)

title('Final Chord Distribution');
xlabel('c/R');
ylabel('r/R');
hold on
plot(r_ad_lin,cv_ad_lin,'r','LineWidth',2)
legend('Optimal Chord','Corrected and Linearized Chord')

%%

cv_final = cv_ad_lin.*R;
r_ad_final = r_ad_tot;

%% a, ap correction with chord known
a_corrected=[];
ap_corrected=[];
ac = a_gv(end);


%% a and ap corrected w/ ONLY A PROFILE EFFECT
r_ad_A = [rh:drad:ra]';
e=0.5;

%while err > tol
    x = rh*lam;
    [kh] = find_tol (xv,x,tol_h);
    
    a = a_gv(kh);
    ap= ap_gv(kh);

  err1= tol +1;
  err2= tol +1;

  while err1 > tol || err2> tol

    x = rh*lam;
    %fi = atan(ap*x/a);
    fi = atan((1-a)/((1+ap)*x));
    w_inf = V0*(1-a)/sin(fi);
    
    c=cv_final(1);
    Re=c*w_inf/v;
    
    [av,lv,dv,Ev,Emax] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
    
    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
    b_c=rad2deg(fi)-aopt;
    F = (2/pi)*acos(exp(Nb*(rh-1)/(2*rh*sin(fi))));
    
    sigma=c*Nb/(2*pi*rh*R);

    if a>ac
        
        func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(lopt*cos(fi)+dopt*sin(fi))*sigma*w_inf^2);
       
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(lopt*sin(fi)-dopt*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);

    else
        a_n =sigma*w_inf*(lopt*cos(fi)+dopt*sin(fi))/(4*V0*F*sin(fi));
        ap_n=sigma*w_inf*(lopt*sin(fi)-dopt*cos(fi))/(4*V0*F*x*sin(fi));
    end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n; 
    ap=(1-e)*ap+(e)*ap_n;
    
  end

a_corrected(1) =a;
ap_corrected(1)=ap;
Fv(1)=F;

aopt_v(1) = aopt;
lopt_v(1) = lopt;
dopt_v(1) = dopt;
w_inf_v(1) = w_inf;
Re_v(1) = Re;
bc_v(1) = b_c;
fi_v_corrected(1)=fi;

for j = 2:length(r_ad_A)
    x = r_ad_A(j)*lam;

    [kh] = find_tol (xv,x,tol_h);
    
    a = a_gv(kh);
    ap= ap_gv(kh);

    err1= tol +1;
    err2= tol +1;
    
    while err1 > tol || err2> tol
    
    x = r_ad_A(j)*lam;
    %fi = atan(ap*x/a);     
    fi = atan((1-a)/((1+ap)*x));   
    w_inf = V0*(1-a)/sin(fi);
    
    [kh] = find_tol (r_ad_tot,r_ad_A(j),tol_a);

    c=cv_final(kh);
    Re=c*w_inf/v;
        
    [av_A,lv_A,dv_A,Ev_A,Emax_A] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
  
    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);

    b_c=rad2deg(fi)-aopt;   
    F = (2/pi)*acos(exp(Nb*(r_ad_A(j)-1)/(2*r_ad_A(j)*sin(fi))));


    sigma=c*Nb/(2*pi*r_ad_A(j)*R);

     if a>ac
        
         func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(lopt*cos(fi)+dopt*sin(fi))*sigma*w_inf^2);
        
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(lopt*sin(fi)-dopt*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
    
     else
         a_n=sigma*w_inf*(lopt*cos(fi)+dopt*sin(fi))/(4*V0*F*sin(fi));
         ap_n=sigma*w_inf*(lopt*sin(fi)-dopt*cos(fi))/(4*V0*F*x*sin(fi));
     end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;

    end

a_corrected(j)=a;
ap_corrected(j)=ap;
Fv(j)=F;
bc_v(j) = b_c;
aopt_v(j) = aopt;
lopt_v(j) = lopt;
dopt_v(j) = dopt;
w_inf_v(j) = w_inf;

fi_v_corrected(j)=fi; 
Re_v(j) = Re;
end

ja=length(a_corrected);

%% a and ap corrected w/ EFFECT OF A & B

r_ad_AB = [ra+drad:drad:rb]';

for j = 1:length(r_ad_AB)
    
    x = r_ad_AB(j)*lam;
    [kh] = find_tol (xv,x,tol_h);
    
    a = a_gv(kh);
    ap= ap_gv(kh);

    err1= tol +1;
    err2= tol +1;
    
    while err1 > tol || err2> tol
    
    x = r_ad_AB(j)*lam;
    %fi = atan(ap*x/a);     
    fi = atan((1-a)/((1+ap)*x));    
    w_inf = V0*(1-a)/sin(fi);
    
    [kh] = find_tol (r_ad_tot,r_ad_AB(j),tol_a);
    
    c=cv_final(kh);
    Re=c*w_inf/v;
         
    [av_A,lv_A,dv_A,Ev_A,Emax_A] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
    [av_B,lv_B,dv_B,Ev_B,Emax_B] = B_PRO(Re,B_200K,B_500K,B_1000K,B_2000K,B_5000K);
    
    [av, lv, dv, Ev,Emax] = rad_interp(av_A,lv_A,dv_A,Ev_A,av_B,lv_B,dv_B,Ev_B,r_ad_AB(j),ra,rb);
   
    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
    b_c=rad2deg(fi)-aopt;
    F = (2/pi)*acos(exp(Nb*(r_ad_AB(j)-1)/(2*r_ad_AB(j)*sin(fi))));
    
    sigma=c*Nb/(2*pi*r_ad_AB(j)*R);
   
    if a>ac

        func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(lopt*cos(fi)+dopt*sin(fi))*sigma*w_inf^2);
       
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(lopt*sin(fi)-dopt*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
    
    else
        a_n =sigma*w_inf*(lopt*cos(fi)+dopt*sin(fi))/(4*V0*F*sin(fi));
        ap_n=sigma*w_inf*(lopt*sin(fi)-dopt*cos(fi))/(4*V0*F*x*sin(fi));
    end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
    
    end

a_corrected(j+ja)=a;
ap_corrected(j+ja)=ap;
Fv(j+ja)=F; 
bc_v(j+ja) = b_c;
aopt_v(j+ja) = aopt;
lopt_v(j+ja) = lopt;
dopt_v(j+ja) = dopt;
w_inf_v(j+ja) = w_inf;
Re_v(j+ja) = Re;
fi_v_corrected(j+ja)=fi;

end

jab=length(ap_corrected);

%% a and ap corrected w/ EFFECT OF B & C

r_ad_BC = [rb+drad:drad:rc]';

for j = 1:length(r_ad_BC)
    
    x = r_ad_BC(j)*lam;
    [kh] = find_tol (xv,x,tol_h);
    
    a = a_gv(kh);
    ap= ap_gv(kh);
    err1= tol +1;
    err2= tol +1;
    
    while err1 > tol || err2> tol
    
    x = r_ad_BC(j)*lam;
    %fi = atan(ap*x/a);     
    fi = atan((1-a)/((1+ap)*x));
    w_inf = V0*(1-a)/sin(fi);
    
    [kh] = find_tol (r_ad_tot,r_ad_BC(j),tol_a);
    
    c=cv_final(kh);
    Re=c*w_inf/v;
    
    [av_B,lv_B,dv_B,Ev_B,Emax_B] = B_PRO(Re,B_200K,B_500K,B_1000K,B_2000K,B_5000K);
    [av_C,lv_C,dv_C,Ev_C,Emax_C] = C_PRO(Re,C_50K,C_100K,C_200K,C_500K,C_1000K);
    
    [av, lv, dv, Ev,Emax] = rad_interp(av_B,lv_B,dv_B,Ev_B,av_C,lv_C,dv_C,Ev_C,r_ad_BC(j),rb,rc);
    
    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
    b_c=rad2deg(fi)-aopt;

    F = (2/pi)*acos(exp(Nb*(r_ad_BC(j)-1)/(2*r_ad_BC(j)*sin(fi))));

    sigma=c*Nb/(2*pi*r_ad_BC(j)*R);
    
    if a>ac

        func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(lopt*cos(fi)+dopt*sin(fi))*sigma*w_inf^2);
        
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(lopt*sin(fi)-dopt*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
   
    else
        a_n=sigma*w_inf*(lopt*cos(fi)+dopt*sin(fi))/(4*V0*F*sin(fi));
        ap_n=sigma*w_inf*(lopt*sin(fi)-dopt*cos(fi))/(4*V0*F*x*sin(fi));
    end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
   
    end

a_corrected(j+jab)=a;
ap_corrected(j+jab)=ap;
Fv(j+jab)=F; 
bc_v(j+jab) = b_c;
aopt_v(j+jab) = aopt;
lopt_v(j+jab) = lopt;
dopt_v(j+jab) = dopt;
w_inf_v(j+jab) = w_inf;
Re_v(j+jab) = Re;

fi_v_corrected(j+jab)=fi;
end

jbc=length(ap_corrected);

%%

r_ad_C=[rc+drad:drad:rt]';

 for j = 1:length(r_ad_C)
    x = r_ad_C(j)*lam;
    [kh] = find_tol (xv,x,tol_h);
    
    a = a_gv(kh);
    ap= ap_gv(kh);
    err1= tol +1;
    err2= tol +1;
    
    while err1 > tol || err2> tol
    
    x = r_ad_C(j)*lam;
    %fi = atan(ap*x/a);     
    fi = atan((1-a)/((1+ap)*x));    
    w_inf = V0*(1-a)/sin(fi);
    
    [kh] = find_tol (r_ad_tot,r_ad_C(j),tol_a);
    
    c=cv_final(kh);
    Re=c*w_inf/v;
    
    [av_C,lv_C,dv_C,Ev_C,Emax_C] = C_PRO(Re,C_50K,C_100K,C_200K,C_500K,C_1000K);
    
    [aopt, lopt, dopt] = finder(av_C, lv_C, dv_C, Ev_C, Emax_C);
    b_c=rad2deg(fi)-aopt;  

    F = (2/pi)*acos(exp(Nb*(r_ad_C(j)-1)/(2*r_ad_C(j)*sin(fi))));

    sigma=c*Nb/(2*pi*r_ad_C(j)*R);
    
    if a>ac
        
        func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(lopt*cos(fi)+dopt*sin(fi))*sigma*w_inf^2);
    
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(lopt*sin(fi)-dopt*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
    
    else 
        a_n =sigma*w_inf*(lopt*cos(fi)+dopt*sin(fi))/(4*V0*F*sin(fi));
        ap_n=sigma*w_inf*(lopt*sin(fi)-dopt*cos(fi))/(4*V0*F*x*sin(fi));
    end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;

    end

a_corrected(j+jbc)=a;
ap_corrected(j+jbc)=ap;
Fv(j+jbc)=F;
bc_v(j+jbc) = b_c;
aopt_v(j+jbc) = aopt;
lopt_v(j+jbc) = lopt;
dopt_v(j+jbc) = dopt;
w_inf_v(j+jbc) = w_inf;
Re_v(j+jbc) = Re;
fi_v_corrected(j+jbc)=fi;   
 end
%%

 for i=1:length(a_corrected)
     a=a_corrected(i);
        if a<ac
            Cx(i)=4*a*(1-a);
        else
             Cx(i)=4*(a*(1-a)+(a-ac)/4*((a-ac)^2/ac+2*(a-ac)+ac));
        end
 end
 figure;
 Cx_ve();
 hold on;
 plot(a_corrected,Cx,LineWidth=1.5);
 vert=0:0.001:2;
 xx=linspace(min(a_corrected),min(a_corrected),length(vert));
 hold on;
 plot(xx,vert,'m',LineWidth=1.5);
 xx=linspace(max(a_corrected),max(a_corrected),length(vert));
 hold on;
 plot(xx,vert,'m',LineWidth=1.5);

%r_ad=[rh:drad:rt];

%% cp calcs
r_ad_tot = r_ad_tot';
int1=Fv.*ap_corrected.*(1-a_corrected).*(r_ad_tot*lam).^3;
cp1=trapz(r_ad_tot*lam,int1)*8/lam^2;

%%
cv_final = cv_final';
%int2 = (lopt_v.*sin(fi_v_corrected)-dopt_v.*cos(fi_v_corrected)).*w_inf_v.^2/(V0^3)*(3.*cv_final.*r_ad_tot*R*rot_s)/(pi*R^2);
int2 = (lopt_v.*sin(fi_v_corrected)-dopt_v.*cos(fi_v_corrected)).*(w_inf_v.^2)*3.*cv_final.*r_ad_tot*R*rot_s;
cp2 = trapz(r_ad_tot*R,int2)/(V0^3)/(pi*R^2);

%%

Mh=86;
cp2 = gauss_quad(0, lam, int1, Mh)*8/lam^2;

%% FIGURE

figure 
plot(r_ad_tot,a_corrected,'k','LineWidth',2)
hold on

plot(r_ad_tot,ap_corrected,'b','LineWidth',2)
%hold on
%plot(r_ad_tot,linspace(1/3,1/3,length(r_ad_tot)),'LineWidth',2)
%plot(r_ad_tot,linspace(0.3328,0.3328,length(r_ad_tot)))
title('a and ap corrected')

xlabel('r/R');
ylabel('a and ap');
hold on;
yyaxis left;
plot(xv,a_gv,'LineWidth',2);
%ylabel('a[-]');
hold on;
yyaxis right;
plot(xv,ap_gv,'LineWidth',2);
%xlabel('x[-]');
%ylabel('ap[-]');
legend('a corrected', 'ap_corrected','a','a''');
%title('Optimal Glauert Distribution');

%%
figure
plot(r_ad_tot, bc_v,'LineWidth',2)
hold on;
plot(r_ad_tot, aopt_v,'LineWidth',2)
hold on;
plot(r_ad_tot, fi_v_corrected.*180/pi(),'LineWidth',2)
hold on;
legend('βc', 'βinf', 'phi' )
xlabel('r/R');
ylabel('phi βc Beta_inf');
title('Flow, pitch and AoA varition along the blade')

figure
yyaxis left;
plot(xv,a_gv,'LineWidth',2);
%ylabel('a[-]');
hold on;
yyaxis right;
plot(xv,ap_gv,'LineWidth',2);
%xlabel('x[-]');
%ylabel('ap[-]');
legend('a','a''');
title('Optimal Glauert Distribution');


figure
plot(r_ad_tot, Re_v,'LineWidth',2)
legend('Re' )
xlabel('r/R');
ylabel('Re');
title('Reynold Number Varition along the blade')


%%
figure
plot(r_ad_tot, w_inf_v./V0,'LineWidth',2)
xlabel('r/R');
ylabel('Winf/V0');
title('Velocity Ratio along the blade')


%% POWER AND MOMENTUM

P_av = 0.5 * ro * V0^3 * pi * R^2;
cp = cp1;
P = cp*P_av;
cm = cp / lam;

aopt_v_deg=aopt_v.*180/pi;


