% PCA: Análise de Componentes Principais para fitting linear

% Intervalo de interesse
t = 0:pi/200:pi/2;

% Variáveis
x = cos(t);
y = sin(t);

% Criar a matriz A de dados: observações em colunas
A = [x; y];

% Centralizar dados
A = A - mean(A, 2); % média ao longo das linhas

% Matriz de covariância: Cx = AA^T/n
C = A*A'/max(size(x));

% Calcular o autovalor dominante da matriz de covariãncia
% Como queremos apenas um fitting linear, usamos apenas o autovetor
% associado ao primeiro autovalor

tol = 0.000001;

% Método das Potências
[eigenvalue, eigenvector, k, err_pot] = potencias(C, tol);
v_potencias = eigenvector;

% Método de Francis
[autovetores, autovalores, err_fr] = francis(C, tol);
v_francis = autovetores(:, 1); % primeiro autovetor

% PLOT
hold on

% Plot dos dados centralizados
plot(A(1, :), A(2, :));

% Variáveis do fitting linear
t = t - mean(t); % centralizando intervalo de interesse 
x_potencias = t * v_potencias(1);
y_potencias = t * v_potencias(2);
x_francis = t * v_francis(1);
y_francis = t * v_francis(2);

% Plot do fitting linear
plot(x_potencias, y_potencias, 'r');
plot(x_francis, y_francis, 'b--');

legend('Dados', 'Fitting Potências', 'Fitting Francis');
xlabel('Eixo X');
ylabel('Eixo Y');
xlim([-0.75, 0.75]);
ylim([-0.75, 0.75]);
hold off

% Menor M. Potências
 [lambda,y,k] = potencia_inv(C,tol, 0);


% FUNÇÕES ----------------------------------------------------------
% Método das Potências
function [lambda,y,k, erro] = potencias(A,tol)
    k = 0; kmax = 1000; erro = inf;
    n = size(A,1); y0 = zeros(n,1); y0(1) = 1;
    while (erro>tol && k<kmax)
        x = A*y0;
        y = x/norm(x);
        erro = abs(abs(y0'*y)-1);
        y0 = y; k = k+1;
    end
    lambda = y'*A*y;
end

% Método de Francis
function [V,D, erro] = francis(A,tol)
    n = size(A,1);
    V = eye(n);
    erro = inf;
    while erro>tol
        [Q,R] = qr(A);
        A = R*Q;
        V = V*Q;
        erro = max(max(abs(tril(A,-1))));
    end
    D = diag(A);
end

% M. Potência invertido
function [lambda,y,k] = potencia_inv(A,tol,alpha)
    if(nargin==2) alpha = 0; end
    k = 0; kmax = 1000; erro = inf;
    n = size(A,1); y0 = zeros(n,1); y0(1) = 1;
    B = A - alpha*eye(n);
    [L,U] = lu(B);
    while (erro>tol && k<kmax)
        x = U\(L\y0);
        y = x/norm(x);
        erro = abs(abs(y0'*y)-1);
        y0 = y; k = k+1;
    end
    lambda = y'*A*y;
end