%% Metodo dos gradientes

% Sistema linear: Ax = b
A = [20.8592    0.6954   19.6649  -87.6414;
    0.6954    3.3732   13.2975  -33.6684;
   19.6649   13.2975  252.3526 -339.8287;
  -87.6414  -33.6684 -339.8287  834.4150];
b = [1 2 3 4]';

% Parametros
x0 = [0 0 0 0]';
tol = 1e-10;

% Solucao do sistema linear:
[x, k] = grad(A, b, x0, tol);
disp(x);
disp(k);

function [x, k] = grad(A, b, x0, tol)
% A: matriz SPD
kmax = 10000;
for k=1:kmax
    r = b - A*x0;
    alpha = dot(r, r)/dot(r, A*r);
    x_ant = x0;
    x0 = x0 + alpha*r; 

    % calcular erro relativo
    errorel = norm(x0 - x_ant)/norm(x0);
    if errorel < tol
        x = x0;
        k = k - 1;
        return;
    end
    x = x0; 
end
disp('Erro: o metodo nao converge.');
end