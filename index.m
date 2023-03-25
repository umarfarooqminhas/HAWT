function [k] = index (xv,x,tol)

% this function finds the index k of the vector xv such that
% x = close to xv(k)

% xv = vector to be scanned
% x = value to be found in xv
% k s.t. xv(k) = k 

err = tol + 1;

for j = 1 : length(xv)

    diff =abs(xv(j)-x);

    if diff <= err
        k=j;
        err=diff;
    end

end

%maybe we  can use a function similar to the finder that works with Glauert




