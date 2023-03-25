function [data] = palas_off(data,s, a_first, a_end)
%FROM txt files to Matlab = matrix [alpha; Cl; Cd]
% INPUT data = data loaded from .txt
% step of alpha from data.txt = 0.25 ° 
% s = step  for aplha

%loading
alpha = data(:,1); % alpha vector with step of 0.25
Cl = data (:,2);   % Cl associated for each alpha
Cd = data (:,3);    % Cd associated for each alpha

a = a_first: s : a_end; % new ROW of alpha with step of s


l = spline(alpha,Cl,a); 
d = spline(alpha,Cd,a);

delta_cl=abs(l(end)-l(end-1));
delta_cd=abs(d(end)-d(end-1));
alpha=-5.5,
    while alpha>-60
        l=[l,l(end)-delta_cl];
        d=[d,d(end)+delta_cd];
        a=[a,alpha];
        alpha=alpha-0.5;
    end

data = [a;l;d]';

end