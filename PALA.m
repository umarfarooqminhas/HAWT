function [data] = PALA(data)
%FROM txt files to Matlab = matrix [alpha; Cl; Cd
%  for each profile (data) => alpha = 1 col, Cl = 2 col, Cd = 3 col
alpha = data(:,1);
Cl = data (:,2);
Cd = data (:,3);

s = 0.05; %step %may be an input data

a = alpha; % alpha vector with step of 0.25
l = Cl;    % Cl associated for each alpha
d = Cd;    % Cd associated for each alpha

da = [a(1): s : a(end)]; % new ROW of alpha with step of 0.2

l = spline(a,l,da); %l = spline(a,l,s);
d = spline(a,d,da);

data = [da;l;d]';

end