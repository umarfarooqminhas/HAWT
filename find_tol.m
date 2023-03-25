function [k] = find_tol (xv,x,tol)

% this function finds the index k of the vector xv such that
% x = close to xv(k)

% xv = vector to be scanned
% x = value to be found in xv
% k s.t. xv(k) = k 

k=1;

while abs(xv(k)-x) / tol > 1
    k = k+1;
end

%abs((xv(k)-x)/x) > tol

 