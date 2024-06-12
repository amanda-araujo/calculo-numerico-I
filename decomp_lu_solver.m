
% Solução do sistema linear: ax = b

% Leitura de dados 
%prompt = "Qual a matriz A do sistema linear?";
%A = input(prompt);
%prompt = "Qual o vetor b do sistema linear?";
%b = input(prompt);

A = [31 05 19 99; 12 06 19 40; 05 06 19 64; 13 01 20 00];
b = [10 7 19 98];


% Decomposição LU: A = LU
[L, U] = lu_decomp(A);
% Substituição forward: Ly = b
y = sub_forward(L, b);
% Substituição backward: Ux = y
x = sub_backward(U, y);

% Solução
disp(x);

% FUNÇÕES 
% Etapa substituição progressiva: matriz triangular inferior
function x=sub_forward(L,b)
% L: matriz triangular inferior
% b: termo independente
% x: vetor solucao

n=length(b);
x=zeros(n,1);

for i=1:n
    x(i)=(b(i)-L(i,1:i-1)*x(1:i-1))/L(i,i);
end
end

% Etapa substituição regressiva: matriz triangular superior
function x=sub_backward(U,y)
% U: matriz triangular superior
% y: termo independente
% x: vetor solucao
n=length(y);
x=zeros(n,1);
for i=n:-1:1
    x(i)=(y(i)-U(i,i+1:n)*x(i+1:n))/U(i,i);
end
end

% decomposição LU
function [L,U]=lu_decomp(A)
% A: matriz quadrada
% L, U: matrizes triang. inf. e sup., respectivamente.
n=size(A,1); L=eye(n); U=zeros(n);
for k=1:n
    for j=k:n
        U(k,j)=A(k,j);
        for s=1:k-1
            U(k,j)=U(k,j)-L(k,s)*U(s,j);
        end
    end
    for i=k+1:n
        L(i,k)=A(i,k);
        for s=1:k-1
            L(i,k)=L(i,k)-L(i,s)*U(s,k);
        end
        L(i,k)=L(i,k)/U(k,k);
    end
end
end
