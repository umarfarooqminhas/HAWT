function [cp, P_CD] = off_CD (bc_on, D_bc , lam, V0, a_corrected, ap_corrected, cv_final)

% DATA & assumptions

s = 0.5;
a_first = -5;
a_end = 15;

%ROOT

load A_sp4721_200K.txt
load A_sp4721_500K.txt
load A_sp4721_1000K.txt
load A_sp4721_2000K.txt
load A_sp4721_5000K.txt

[A_200K_off]  = palas(A_sp4721_200K,s, a_first, a_end);  %A_sp4721_200K_off
[A_500K_off]  = palas(A_sp4721_500K,s, a_first, a_end);  %A_sp4721_500K_off
[A_1000K_off] = palas(A_sp4721_1000K,s, a_first, a_end); %A_sp4721_1000K_off
[A_2000K_off] = palas(A_sp4721_2000K,s, a_first, a_end); %A_sp4721_2000K_off
[A_5000K_off] = palas(A_sp4721_5000K,s, a_first, a_end); %A_sp4721_5000K_off

%PRIMARY

load B_sp4621_200K.txt
load B_sp4621_500K.txt
load B_sp4621_1000K.txt
load B_sp4621_2000K.txt
load B_sp4621_5000K.txt

[B_200K_off]  = palas(B_sp4621_200K,s, a_first, a_end);  %B_sp4621_200K_off
[B_500K_off]  = palas(B_sp4621_500K,s, a_first, a_end);  %B_sp4621_500K_off
[B_1000K_off] = palas(B_sp4621_1000K,s, a_first, a_end); %B_sp4621_1000K_off
[B_2000K_off] = palas(B_sp4621_2000K,s, a_first, a_end); %B_sp4621_2000K_off
[B_5000K_off] = palas(B_sp4621_5000K,s, a_first, a_end); %B_sp4621_5000K_off

%TIP

load C_S828_50K.txt
load C_S828_100K.txt
load C_S828_200K.txt
load C_S828_500K.txt
load C_S828_1000K.txt

[C_50K_off]   = palas(C_S828_50K,s, a_first, a_end);   %C_S828_50K_off
[C_100K_off]  = palas(C_S828_100K,s, a_first, a_end);  %C_S828_100K_off
[C_200K_off]  = palas(C_S828_200K,s, a_first, a_end);  %C_S828_200K_off
[C_500K_off]  = palas(C_S828_500K,s, a_first, a_end);  %C_S828_500K_off
[C_1000K_off] = palas(C_S828_1000K,s, a_first, a_end); %C_S828_1000K_off


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

e=0.5;

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

p = 0.5;

while p > 0.2
a_cor_off=[];
ap_cor_off=[];
Fv_off=[];
fi_v_off=[];
bc_on = bc_on+D_bc;
% a and ap corrected w/ ONLY A PROFILE EFFECT
r_ad_A = [rh:drad:ra]';

cou = 0; %counter

    a = a_corrected(1);
    ap= ap_corrected(1);

  err1= tol +1;
  err2= tol +1;
  err3 = max(err1,err2);
  it = 1;
  it_max = 1000;
  a_vec_it = a;
    
 while err3 > tol && it < it_max
        
        it = it+1;

    x = rh*lam;
    %fi = atan(ap*x/a);     
    fi = atan((1-a)/((1+ap)*x));
    w_inf = V0*(1-a)/sin(fi);

    b_c = bc_on(1);
    a_off = rad2deg(fi)-b_c;
    
    c=cv_final(1);
    Re=c*w_inf/v;
     
    AR = pi*R/c;
    if AR<=50
        CD_max=1.11+0.018*AR;
    else
        CD_max=2.01;
    end

    [av,lv,dv,Ev,Emax] = A_PRO(Re,A_200K_off,A_500K_off,A_1000K_off,A_2000K_off,A_5000K_off);
        
    if a_off>a_first && a_off<a_end
        [k] = find_tol (av,a_off,tol_alp);
        l_off = lv(k);
        d_off = dv(k);
        E_off = l_off/d_off;
    else
        [k] = find_tol(lv,max(lv),0.01);
        a_stall=av(k);
        cl_stall=lv(k);
        cd_stall=dv(k);
    
        K_L=(cl_stall-CD_max*sind(a_stall)*cosd(a_stall))*sind(a_stall)/(cosd(a_stall))^2;
        K_D=(cd_stall-CD_max*(sind(a_stall))^2)/cosd(a_stall);
    
        l_off=CD_max/2*sind(2*a_off)+K_L*(cosd(a_off))^2/sind(a_off);
        d_off=CD_max*(sind(a_off))^2+K_D*cosd(a_off);
        l_off_v(1)=l_off;
        d_off_v(1)=d_off;
    end

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
    err3 = max(err1,err2);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
    a_vec_it = [a_vec_it; a];
    
  end

  [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
  
  if E_off < 0.95*Emax && a_off > aopt
      cou = cou +1;
  end

a1 = a_vec_it(end);
a2 = a_vec_it(end-1);
a = (a1+a2)/2;
a_cor_off(1) =a;
ap_cor_off(1)=ap;
Fv_off(1)=F;
fi_v_off(1)=fi;


for j = 2:length(r_ad_A)
    
    a = a_corrected(j);
    ap= ap_corrected(j);

    err1= tol +1;
    err2= tol +1;
    err3 = max(err1,err2);
    it = 1;
    it_max = 1000;
    a_vec_it = a;
    
 while err3 > tol && it < it_max
        
        it = it+1;
    x = r_ad_A(j)*lam;
    %fi = atan(ap*x/a);    
    fi = atan((1-a)/((1+ap)*x));
    w_inf = V0*(1-a)/sin(fi);

    b_c = bc_on(j);
    a_off = rad2deg(fi)-b_c;
    
    c=cv_final(j);
    Re=c*w_inf/v;
     
    AR = pi*R/c;
    if AR<=50
        CD_max=1.11+0.018*AR;
    else
        CD_max=2.01;
    end

    [av,lv,dv,Ev,Emax] = A_PRO(Re,A_200K_off,A_500K_off,A_1000K_off,A_2000K_off,A_5000K_off);
              
    if a_off>a_first && a_off<a_end    
        [k] = find_tol (av,a_off,tol_alp);
        l_off = lv(k);
        d_off = dv(k);
        E_off = l_off/d_off;
    else
        [k] = find_tol(lv,max(lv),0.01);
        a_stall=av(k);
        cl_stall=lv(k);
        cd_stall=dv(k);
    
        K_L=(cl_stall-CD_max*sind(a_stall)*cosd(a_stall))*sind(a_stall)/(cosd(a_stall))^2;
        K_D=(cd_stall-CD_max*(sind(a_stall))^2)/cosd(a_stall);
    
        l_off=CD_max/2*sind(2*a_off)+K_L*(cosd(a_off))^2/sind(a_off);
        d_off=CD_max*(sind(a_off))^2+K_D*cosd(a_off);
        l_off_v(j)=l_off;
        d_off_v(j)=d_off;
    end   
    
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
    err3 = max(err1,err2);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
    a_vec_it = [a_vec_it; a];

    end

    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
  
  if E_off < 0.95*Emax && a_off > aopt
      cou = cou +1;
  end

a1 = a_vec_it(end);
a2 = a_vec_it(end-1);
a = (a1+a2)/2;
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
    err3 = max(err1,err2);
    it = 1;
    it_max = 1000;
    a_vec_it = a;
    
 while err3 > tol && it < it_max
        
        it = it+1;
    x = r_ad_AB(j)*lam;
    %fi = atan(ap*x/a);    
    fi = atan((1-a)/((1+ap)*x));
    w_inf = V0*(1-a)/sin(fi);

    b_c = bc_on(j+ja);
    a_off = rad2deg(fi)-b_c;
    c=cv_final(j+ja);
    Re=c*w_inf/v;
     
    AR = pi*R/c;
    if AR<=50
        CD_max=1.11+0.018*AR;
    else
        CD_max=2.01;
    end

    [av_A,lv_A,dv_A,Ev_A,Emax_A] = A_PRO(Re,A_200K_off,A_500K_off,A_1000K_off,A_2000K_off,A_5000K_off);
    [av_B,lv_B,dv_B,Ev_B,Emax_B] = B_PRO(Re,B_200K_off,B_500K_off,B_1000K_off,B_2000K_off,B_5000K_off);
        
    [av, lv, dv, Ev,Emax] = rad_interp(av_A,lv_A,dv_A,Ev_A,av_B,lv_B,dv_B,Ev_B,r_ad_AB(j),ra,rb);
               
    if a_off>a_first && a_off<a_end    
        [k] = find_tol (av,a_off,tol_alp);
        l_off = lv(k);
        d_off = dv(k);    
        E_off = l_off/d_off;
    else
        [k] = find_tol(lv,max(lv),0.01);
        a_stall=av(k);
        cl_stall=lv(k);
        cd_stall=dv(k);
    
        K_L=(cl_stall-CD_max*sind(a_stall)*cosd(a_stall))*sind(a_stall)/(cosd(a_stall))^2;
        K_D=(cd_stall-CD_max*(sind(a_stall))^2)/cosd(a_stall);
    
        l_off=CD_max/2*sind(2*a_off)+K_L*(cosd(a_off))^2/sind(a_off);
        d_off=CD_max*(sind(a_off))^2+K_D*cosd(a_off);
        l_off_v(j)=l_off;
        d_off_v(j)=d_off;
    end

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
    err3 = max(err1,err2);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
    a_vec_it = [a_vec_it; a];
    
    end

    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
  
  if E_off < 0.95*Emax && a_off > aopt
      cou = cou +1;
  end

a1 = a_vec_it(end);
a2 = a_vec_it(end-1);
a = (a1+a2)/2;
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
    err3 = max(err1,err2);
    it = 1;
    it_max = 1000;
    a_vec_it = a;
    
    while err3 > tol && it < it_max
    
        it = it+1;
    x = r_ad_BC(j)*lam;
    %fi = atan(ap*x/a);    
    fi = atan((1-a)/((1+ap)*x));
    w_inf = V0*(1-a)/sin(fi);

    b_c = bc_on(j+jab);
    a_off = rad2deg(fi)-b_c;
    
    c=cv_final(j+jab);
    Re=c*w_inf/v;
     
    AR = pi*R/c;
    if AR<=50
        CD_max=1.11+0.018*AR;
    else
        CD_max=2.01;
    end

    [av_B,lv_B,dv_B,Ev_B,Emax_B] = B_PRO(Re,B_200K_off,B_500K_off,B_1000K_off,B_2000K_off,B_5000K_off);
    [av_C,lv_C,dv_C,Ev_C,Emax_C] = C_PRO(Re,C_50K_off,C_100K_off,C_200K_off,C_500K_off,C_1000K_off);
        
    [av, lv, dv, Ev,Emax] = rad_interp(av_B,lv_B,dv_B,Ev_B,av_C,lv_C,dv_C,Ev_C,r_ad_BC(j),rb,rc);
                
    if a_off>a_first && a_off<a_end
        [k] = find_tol (av,a_off,tol_alp);
        l_off = lv(k);
        d_off = dv(k);   
        E_off = l_off/d_off;
    else
        [k] = find_tol(lv,max(lv),0.01);
        a_stall=av(k);
        cl_stall=lv(k);
        cd_stall=dv(k);
    
        K_L=(cl_stall-CD_max*sind(a_stall)*cosd(a_stall))*sind(a_stall)/(cosd(a_stall))^2;
        K_D=(cd_stall-CD_max*(sind(a_stall))^2)/cosd(a_stall);
    
        l_off=CD_max/2*sind(2*a_off)+K_L*(cosd(a_off))^2/sind(a_off);
        d_off=CD_max*(sind(a_off))^2+K_D*cosd(a_off);
        l_off_v(j)=l_off;
        d_off_v(j)=d_off;
    end

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
    err3 = max(err1,err2);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
    a_vec_it = [a_vec_it; a];
   
    end

    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
  
  if E_off < 0.95*Emax && a_off > aopt
      cou = cou +1;
  end

a1 = a_vec_it(end);
a2 = a_vec_it(end-1);
a = (a1+a2)/2;
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
    err3 = max(err1,err2);
    it = 1;
    it_max = 1000;
    a_vec_it = a;
    
    while it < it_max && err3 > tol

        it = it + 1;

        x = r_ad_C(j)*lam;
        %fi = atan(ap*x/a);   
        fi = atan((1-a)/((1+ap)*x));
        w_inf = V0*(1-a)/sin(fi);
    
        b_c = bc_on(j+jbc);
        a_off = rad2deg(fi)-b_c;
            
        c=cv_final(j+jbc);
        Re=c*w_inf/v;
         
    AR = pi*R/c;
    if AR<=50
        CD_max=1.11+0.018*AR;
    else
        CD_max=2.01;
    end

    [av,lv,dv,Ev,Emax] = C_PRO(Re,C_50K_off,C_100K_off,C_200K_off,C_500K_off,C_1000K_off);

    if a_off>a_first && a_off<a_end
        [k] = find_tol (av,a_off,tol_alp);
        l_off = lv(k);
        d_off = dv(k);  
        E_off = l_off/d_off;
    else
        [k] = find_tol(lv,max(lv),0.01);
        a_stall=av(k);
        cl_stall=lv(k);
        cd_stall=dv(k);
    
        K_L=(cl_stall-CD_max*sind(a_stall)*cosd(a_stall))*sind(a_stall)/(cosd(a_stall))^2;
        K_D=(cd_stall-CD_max*(sind(a_stall))^2)/cosd(a_stall);
    
        l_off=CD_max/2*sind(2*a_off)+K_L*(cosd(a_off))^2/sind(a_off);
        d_off=CD_max*(sind(a_off))^2+K_D*cosd(a_off);
        l_off_v(j)=l_off;
        d_off_v(j)=d_off;
    end

        F = (2/pi)*acos(exp(Nb*(r_ad_C(j)-1)/(2*r_ad_C(j)*sin(fi))));

        sigma=c*Nb/(2*pi*r_ad_C(j)*R);

        if a>ac
        
          func = @ (an) (4*(an*(1-an)+(an-ac)/4*((an-ac)^2/ac+2*(an-ac)+ac))*V0^2*F-(l_off*cos(fi)+d_off*sin(fi))*sigma*w_inf^2);
          %func_w = @ (an) (4*(an*(1-an)+(an-ac)^2)*V0^2*F-(l_off*cos(fi)+d_off*sin(fi))*sigma*w_inf^2);
    
          a0=a; %(1-0.8)*a+(0.8)*a_n
          a_n=fzero(func,a0);
          ap_n=(l_off*sin(fi)-d_off*cos(fi))*sigma*w_inf^2/(4*(1-a_n)*F*x*V0^2);
    
        else 
          a_n =sigma*w_inf*(l_off*cos(fi)+d_off*sin(fi))/(4*V0*F*sin(fi));
          ap_n=sigma*w_inf*(l_off*sin(fi)-d_off*cos(fi))/(4*V0*F*x*sin(fi));
        end
    
    err1 = abs(a - a_n);
    err2 = abs(ap - ap_n);
    err3 = max(err1,err2);
    a=(1-e)*a+(e)*a_n;
    ap=(1-e)*ap+(e)*ap_n;
    a_vec_it = [a_vec_it; a];        

    end

    a1 = a_vec_it(end);
    a2 = a_vec_it(end-1);
    a = (a1+a2)/2;

    [aopt, lopt, dopt] = finder(av, lv, dv, Ev, Emax);
  
  if E_off < 0.95*Emax && a_off > aopt
      cou = cou +1;
  end

a_cor_off(j+jbc)=a;
ap_cor_off(j+jbc)=ap;
Fv_off(j+jbc)=F;
fi_v_off(j+jbc)=fi;   
 end
D_bc = D_bc+0.5;

p=cou/length(av);
end 

 %stalled blade

% cp calcs
%%
int=Fv_off.*ap_cor_off.*(1-a_cor_off).*(r_ad_tot*lam).^3;
cp=trapz(r_ad_tot*lam,int)*8/lam^2;

P_av = 0.5 * ro * V0^3 * pi * R^2;
P_CD = cp*P_av;