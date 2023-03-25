function I = trapcomp(a, b, N, f)

% N = numero dei sottointervalli
% N+1 = numero di nodi
% a e b = estremi di integrazione

h = (b-a)/N; % ampiezza dei sottointervalli

I = 0; % inizializzo valore integrale

% ciclo sul numero di sottointervalli
for k = 0:N-1
    
    xs = a +  k    * h;    % estremo sx del sottointervallo    
   
    xd = a + (k+1) * h;    % estremo dx del sottointervallo (1 + N = N+1)

    %I_k = 0.5*(xd-xs) * ( f(xs) + f(xd) );
    %I = I + I_k;

    I = I + 0.5*(xd-xs)* (f(xs)+ f(xd));   
end