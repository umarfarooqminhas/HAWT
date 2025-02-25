function [p] = Cx_ve()
at=1/3;
step_a=0.001;

a_v_1 = 0:step_a:at;
a_v_2 = at+step_a:step_a:1;

for j1 = 1:length(a_v_1)
    a = a_v_1(j1);
    cx_1(j1) = 4*(a*(1-a));
end

for j2 = 1:length(a_v_2)
    a = a_v_2(j2);
    cx_2(j2) = 4*(a*(1-a)+(a-at)/4*((a-at)^2/at+2*(a-at)+at));
end

a_v = [a_v_1,a_v_2]';

Cx_v = [cx_1,cx_2]';

p=plot(a_v,Cx_v,LineWidth=1.5);
end