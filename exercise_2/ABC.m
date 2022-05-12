clear all
close all
clc

data = readtable('data_lineal_3.csv'); % 1, 2, 3

X = data.x;
Y = data.y;

n = size(X,1);
X = [ones(n,1) X];

f = @(w) (1/(2*n))*sum((Y-(X*w)).^2);

xl = [-10; -10];
xu = [10; 10];

G = 150;
N = 50;
D = 2;
L = 15; %limite de intentos menor a iteraciones

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

w = x(:,igb);
disp([' Coeficientes: w = [' num2str(w') '], f(w) = ' num2str(f(w))])

Yp = X*w;
R2 = 1 - sum((Y-Yp).^2)/sum((Y-mean(Y)).^2);
disp([' R2 score = ' num2str(R2)])

figure
hold on
grid on
plot(X(:,2),Y,'bo','LineWidth',2,'MarkerSize',10)
plot(X(:,2),Yp,'g-','LineWidth',3,'MarkerSize',10)
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
legend({'Muestras','Regresion'},'FontSize',13)


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