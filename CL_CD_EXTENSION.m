clear all
Re=5*10^6;  
s=0.05;
a_first=-5;
a_end=15;
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


[av,lv,dv,Ev,Emax] = A_PRO(Re,A_200K_off,A_500K_off,A_1000K_off,A_2000K_off,A_5000K_off);
alpha=-60:0.5:-5;
for j=1:length(alpha)
    a_off=alpha(j);
     [k] = find_tol(lv,max(lv),0.01);
     a_stall=av(k);
     cl_stall=lv(k);
     cd_stall=dv(k);
     CD_max=2.01;
     K_L=(cl_stall-CD_max*sind(a_stall)*cosd(a_stall))*sind(a_stall)/(cosd(a_stall))^2;
     K_D=(cd_stall-CD_max*(sind(a_stall))^2)/cosd(a_stall);
     l_off=CD_max/2*sind(2*a_off)+K_L*(cosd(a_off))^2/sind(a_off);
     
     d_off=CD_max*(sind(a_off))^2+K_D*cosd(a_off);
     l_off_v(j)=l_off;
     d_off_v(j)=d_off;
end
%plot(alpha,l_off_v)
%hold on;
%plot(alpha,d_off_v)
%hold on;
plot(alpha,l_off_v./d_off_v)
%legend('CL','CD')
