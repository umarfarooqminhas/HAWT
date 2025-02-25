function [cp, P] = off_AC (bc_on, lam, V0, a_corrected, ap_corrected, cv_final)

% DATA & assumptions

%ROOT
s=0.05;

load A_sp4721_200K.txt;
load A_sp4721_500K.txt;
load A_sp4721_1000K.txt;
load A_sp4721_2000K.txt;
load A_sp4721_5000K.txt;

[A_200K]  = palas(A_sp4721_200K,s);  %A_sp4721_200K
[A_500K]  = palas(A_sp4721_500K,s);  %A_sp4721_500K
[A_1000K] = palas(A_sp4721_1000K,s); %A_sp4721_1000K
[A_2000K] = palas(A_sp4721_2000K,s); %A_sp4721_2000K
[A_5000K] = palas(A_sp4721_5000K,s); %A_sp4721_5000K

%PRIMARY
s=0.05;

load B_sp4621_200K.txt;
load B_sp4621_500K.txt;
load B_sp4621_1000K.txt;
load B_sp4621_2000K.txt;
load B_sp4621_5000K.txt;

[B_200K]  = palas(B_sp4621_200K,s);  %B_sp4621_200K
[B_500K]  = palas(B_sp4621_500K,s);  %B_sp4621_500K
[B_1000K] = palas(B_sp4621_1000K,s); %B_sp4621_1000K
[B_2000K] = palas(B_sp4621_2000K,s); %B_sp4621_2000K
[B_5000K] = palas(B_sp4621_5000K,s); %B_sp4621_5000K

%TIP
s=0.05;

load C_S828_50K.txt;
load C_S828_100K.txt;
load C_S828_200K.txt;
load C_S828_500K.txt;
load C_S828_1000K.txt;

[C_50K]   = palas(C_S828_50K,s);   %C_S828_50K
[C_100K]  = palas(C_S828_100K,s);  %C_S828_100K
[C_200K]  = palas(C_S828_200K,s);  %C_S828_200K
[C_500K]  = palas(C_S828_500K,s);  %C_S828_500K
[C_1000K] = palas(C_S828_1000K,s); %C_S828_1000K


D=110;                    % DIAMETER	                [m]
R = D/2;                  % Radius                      [m]
Nb=3;	                  % Number of Blades            [-]		
mi=1.81*10^(-5);          % Dynamic viscosity           [Kg/ms]
ro=1.22;                  % Rho = density               [kg/m3]
v=mi / ro;                % Kinematic Viscosity         [m2/s]
z = 100;                  % altidude	                [m]
%V0 = 7.30;                % Wind velocity (undisturbed) [m/s]
%lam = 6.90;               % lambda = tip speed ratio    [-]

% tolerances selection

drad = 0.01;          % 0.01 0.001

tol = 0.01;           % 0.01

tol_alp = 0.4; % tolerance to find the alpha off design

rh = 0.10; % adimensional radius ratio - hub
ra = 0.20; % ROOT PROFILE = A
rr = 0.40; % adimensional radius ratio - root
rb = 0.50; % PRIMARY PROFILE = B
rp = 0.80; % adimensional radius ratio - primary
rc = 0.86; % TIP PROFILE = C
rt = 0.95; % adimensional radius ratio - tip

r_ad_tot = [rh:drad:rt];

%rf = r_ad(end); % adimensional radius ratio - final

R_h = rh * R ; % R_hub
R_r = rr * R ; % R_root
R_p = rp * R ; % R_primary
R_t = rt * R ; % R_tip
%R_f = rf * R ; % R_final
ac = 0.3328;
% a and ap corrected w/ ONLY A PROFILE EFFECT
r_ad_A = [rh:drad:ra]';
e=0.5;

    a = a_corrected(1);
    ap= ap_corrected(1);

  err1= tol +1;
  err2= tol +1;

  while err1 > tol || err2> tol

    x = rh*lam;
    fi = atan(ap*x/a);
    w_inf = V0*(1-a)/sin(fi);

    b_c = bc_on(1);
    a_off = rad2deg(fi)-b_c;
    
    c=cv_final(1);
    Re=c*w_inf/v;
    
    [av,lv,dv,Ev,Emax] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
    
    [k] = find_tol (av,a_off,tol_alp);
    l_off = lv(k);
    d_off = dv(k);
    
    F = (2/pi)*acos(exp(Nb*(rh-1)/(2*rh*sin(fi))));
    
    sigma=c*Nb/(2*pi*rh*R);

    if a>ac
        
        func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(l_off*cos(fi)+d_off*sin(fi))*sigma*w_inf^2);
       
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(l_off*sin(fi)-d_off*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);

    else
        a_n =sigma*w_inf*(l_off*cos(fi)+d_off*sin(fi))/(4*V0*F*sin(fi));
        ap_n=sigma*w_inf*(l_off*sin(fi)-d_off*cos(fi))/(4*V0*F*x*sin(fi));
    end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n; 
    ap=(1-e)*ap+(e)*ap_n;
    
  end

a_cor_off(1) =a;
ap_cor_off(1)=ap;
Fv_off(1)=F;
fi_v_off(1)=fi;

for j = 2:length(r_ad_A)
    
    a = a_corrected(j);
    ap= ap_corrected(j);

    err1= tol +1;
    err2= tol +1;
    
    while err1 > tol || err2> tol
    
    x = r_ad_A(j)*lam;
    fi = atan(ap*x/a);
    w_inf = V0*(1-a)/sin(fi);

    b_c = bc_on(j);
    a_off = rad2deg(fi)-b_c;
    
    c=cv_final(j);
    Re=c*w_inf/v;
        
    [av_A,lv_A,dv_A,Ev_A,Emax_A] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
  
    [k] = find_tol (av_A,a_off,tol_alp);
    l_off = lv(k);
    d_off = dv(k);
       
    F = (2/pi)*acos(exp(Nb*(r_ad_A(j)-1)/(2*r_ad_A(j)*sin(fi))));

    sigma=c*Nb/(2*pi*r_ad_A(j)*R);

     if a>ac
        
         func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(l_off*cos(fi)+d_off*sin(fi))*sigma*w_inf^2);
        
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(l_off*sin(fi)-d_off*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
    
     else
         a_n=sigma*w_inf*(l_off*cos(fi)+d_off*sin(fi))/(4*V0*F*sin(fi));
         ap_n=sigma*w_inf*(l_off*sin(fi)-d_off*cos(fi))/(4*V0*F*x*sin(fi));
     end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;

    end

a_cor_off(j)=a;
ap_cor_off(j)=ap;
Fv_off(j)=F;
fi_v_off(j)=fi; 

end

ja=length(a_cor_off);

% a and ap corrected w/ EFFECT OF A & B

r_ad_AB = [ra+drad:drad:rb]';

for j = 1:length(r_ad_AB)

    a = a_corrected(j+ja);
    ap= ap_corrected(j+ja);
    
    err1= tol +1;
    err2= tol +1;
    
    while err1 > tol || err2> tol
    
    x = r_ad_AB(j)*lam;
    fi = atan(ap*x/a);
    w_inf = V0*(1-a)/sin(fi);

    b_c = bc_on(j+ja);
    a_off = rad2deg(fi)-b_c;
    
    c=cv_final(j+ja);
    Re=c*w_inf/v;
        
    [av_A,lv_A,dv_A,Ev_A,Emax_A] = A_PRO(Re,A_200K,A_500K,A_1000K,A_2000K,A_5000K);
    [av_B,lv_B,dv_B,Ev_B,Emax_B] = B_PRO(Re,B_200K,B_500K,B_1000K,B_2000K,B_5000K);
    
    [av, lv, dv, Ev,Emax] = rad_interp(av_A,lv_A,dv_A,Ev_A,av_B,lv_B,dv_B,Ev_B,r_ad_AB(j),ra,rb);
   
    [k] = find_tol (av,a_off,tol_alp);
    l_off = lv(k);
    d_off = dv(k);        
    
    F = (2/pi)*acos(exp(Nb*(r_ad_AB(j)-1)/(2*r_ad_AB(j)*sin(fi))));
    
    sigma=c*Nb/(2*pi*r_ad_AB(j)*R);
   
    if a>ac

         func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(l_off*cos(fi)+d_off*sin(fi))*sigma*w_inf^2);
       
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(l_off*sin(fi)-d_off*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
    
    else
        a_n =sigma*w_inf*(l_off*cos(fi)+d_off*sin(fi))/(4*V0*F*sin(fi));
        ap_n=sigma*w_inf*(l_off*sin(fi)-d_off*cos(fi))/(4*V0*F*x*sin(fi));
    end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
    
    end

a_cor_off(j+ja)=a;
ap_cor_off(j+ja)=ap;
Fv_off(j+ja)=F; 
fi_v_off(j+ja)=fi;
end

jab=length(ap_cor_off);

% a and ap corrected w/ EFFECT OF B & C

r_ad_BC = [rb+drad:drad:rc]';

for j = 1:length(r_ad_BC)
    
    a = a_corrected(j+jab);
    ap= ap_corrected(j+jab);

    err1= tol +1;
    err2= tol +1;
    
    while err1 > tol || err2> tol
    
    x = r_ad_BC(j)*lam;
    fi = atan(ap*x/a);
    w_inf = V0*(1-a)/sin(fi);

    b_c = bc_on(j+jab);
    a_off = rad2deg(fi)-b_c;
    
    c=cv_final(j+jab);
    Re=c*w_inf/v;
    
    [av_B,lv_B,dv_B,Ev_B,Emax_B] = B_PRO(Re,B_200K,B_500K,B_1000K,B_2000K,B_5000K);
    [av_C,lv_C,dv_C,Ev_C,Emax_C] = C_PRO(Re,C_50K,C_100K,C_200K,C_500K,C_1000K);
    
    [av, lv, dv, Ev,Emax] = rad_interp(av_B,lv_B,dv_B,Ev_B,av_C,lv_C,dv_C,Ev_C,r_ad_BC(j),rb,rc);
     
    [k] = find_tol (av,a_off,tol_alp);
    l_off = lv(k);
    d_off = dv(k);        
       
    F = (2/pi)*acos(exp(Nb*(r_ad_BC(j)-1)/(2*r_ad_BC(j)*sin(fi))));

    sigma=c*Nb/(2*pi*r_ad_BC(j)*R);
    
    if a>ac

         func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(l_off*cos(fi)+d_off*sin(fi))*sigma*w_inf^2);
        
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(l_off*sin(fi)-d_off*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
   
    else
        a_n=sigma*w_inf*(l_off*cos(fi)+d_off*sin(fi))/(4*V0*F*sin(fi));
        ap_n=sigma*w_inf*(l_off*sin(fi)-d_off*cos(fi))/(4*V0*F*x*sin(fi));
    end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
   
    end

a_cor_off(j+jab)=a;
ap_cor_off(j+jab)=ap;
Fv_off(j+jab)=F; 
fi_v_off(j+jab)=fi;
end

jbc=length(ap_cor_off);

r_ad_C=[rc+drad:drad:rt]';

 for j = 1:length(r_ad_C)
    
    a = a_corrected(j+jbc);
    ap= ap_corrected(j+jbc);
    err1= tol +1;
    err2= tol +1;
    
    while err1 > tol || err2> tol
    
    x = r_ad_C(j)*lam;
    fi = atan(ap*x/a);
    w_inf = V0*(1-a)/sin(fi);
    
    b_c = bc_on(j+jbc);
    a_off = rad2deg(fi)-b_c;
        
    c=cv_final(j+jbc);
    Re=c*w_inf/v;
    
    [av_C,lv_C,dv_C,Ev_C,Emax_C] = C_PRO(Re,C_50K,C_100K,C_200K,C_500K,C_1000K);
    
    [k] = find_tol (av_C,a_off,tol_alp);
    l_off = lv(k);
    d_off = dv(k);  

    F = (2/pi)*acos(exp(Nb*(r_ad_C(j)-1)/(2*r_ad_C(j)*sin(fi))));

    sigma=c*Nb/(2*pi*r_ad_C(j)*R);
    
    if a>ac
        
         func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(l_off*cos(fi)+d_off*sin(fi))*sigma*w_inf^2);
    
        a0=a;
        a_n=fzero(func,a0);
        ap_n=(l_off*sin(fi)-d_off*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
    
    else 
        a_n =sigma*w_inf*(l_off*cos(fi)+d_off*sin(fi))/(4*V0*F*sin(fi));
        ap_n=sigma*w_inf*(l_off*sin(fi)-d_off*cos(fi))/(4*V0*F*x*sin(fi));
    end

    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;

    end

a_cor_off(j+jbc)=a;
ap_cor_off(j+jbc)=ap;
Fv_off(j+jbc)=F;
fi_v_off(j+jbc)=fi;   
 end

% cp calcs
%%
int=Fv_off.*ap_cor_off.*(1-a_cor_off).*(r_ad_tot*lam).^3;
cp=trapz(r_ad_tot*lam,int)*8/lam^2;

P_av = 0.5 * ro * V0^3 * pi * R^2;
P = cp*P_av;