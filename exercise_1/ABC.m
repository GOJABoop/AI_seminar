clear all
close all
clc

%Sistema 1
% A = [2 3; 5 1];
% b = [7; 11];
% D = 2;
% xl = [-10; -10];
% xu = [10; 10];

%Sistema 2
% A = [2 3 4; 1 2 3; 5 1 0; 3 4 1];
% b = [19; 13; 11; 13];
% D = 3;
% xl = [-10; -10; -10];
% xu = [10; 10; 10];

%Sistema 3
% A = [2 3; 5 4; 2 5; 4 1; 0.5 0.5];
% b = [-5; 5; -15; 15; 0];
% D = 2;
% xl = [-10; -10];
% xu = [10; 10];

%Sistema 4
A = [1 -2 2 -3; 3 4 -1 1; 2 -3 2 -1; 1 1 -3 -2];
b = [15; -6; 17; -7];
D = 4;
xl = [-20; -20; -20; -20];
xu = [20; 20; 20; 20];

m = size(b,1);
f = @(x) (1/(2*m))*sum((b-A*x).^2);

G = 200;
N = 50;
L = 40; %limite de intentos menor a iteraciones

Pf = 30;
Po = N-Pf;

x = zeros(D,Pf);
l = zeros(1,Pf);
aptitud = zeros(1,Pf);
fitness = zeros(1,Pf);
f_plot = zeros(1,G);

for i=1:Pf
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = f(x(:,i));
end

for g=1:G 
    %Etapa de abejas empleadas
    for i=1:Pf
        k = i;
        while k == i
            k = randi([1,Pf]);
        end

        j = randi([1,D]);
        phi = 2 * rand()-1;
        v = x(:,i);
        v(j) = x(j,i) + phi*(x(j,i)-x(j,k));
        fv = f(v);

        if fv < fitness(i)
            x(:,i) = v;
            fitness(i) = fv;
            l(i) = 0;
        else
            l(i) = l(i)+1;
        end
        
        if fitness(i) >= 0
            aptitud(i) = 1/(1+fitness(i));
        else
            aptitud(i) = 1 + abs(fitness(i));
        end
    end

    %Etapa de abejas observadoras
    for i=1:Po
        m = Seleccion(aptitud);
        
        k = m;
        while k == m
            k = randi([1,Pf]);
        end
        
        j = randi([1,D]);
        phi = 2 * rand()-1;
   
        v = x(:,m);
        v(j) = x(j,m) + phi *(x(j,m)-x(j,k));
        fv = f(v);

        if fv < fitness(m)
            x(:,m) = v;
            fitness(m) = fv;
            l(m) = 0;
        else
            l(m) = l(m)+1;
        end
    end

    %Etapa de abejas exploradoras
    for i=1:Pf
        if l(i) > L
            x(:,i) = xl + (xu-xl).*rand(D,1);
            fitness(i) = f(x(:,i));
            l(i) = 0;
        end
    end

    f_plot(g) = min(fitness);
end

[~,igb] = min(fitness);

disp(['f(x)= ' num2str(fitness(igb))])
disp('x= ')
disp(x(:,igb))

%validation
disp('Validation')
disp(A*x(:,igb))
disp('b= ')
disp(b)

function[n] = Seleccion(aptitud) %Ruleta
    aptitud_total = sum(aptitud);
    N = numel(aptitud);

    r = rand;
    P_sum = 0;

    for i=1:N
        P_sum = P_sum + aptitud(i)/aptitud_total;
        if P_sum >= r
            n = i;
            return
        end
    end

    n = N;
end