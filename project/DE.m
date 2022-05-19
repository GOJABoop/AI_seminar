clear all
close all
clc

animacion = 0;
pd = [-0.8; 0.4]; %punto

a_1 = 0.35;
a_2 = 0.35;
a_3 = 0.25;

beta = 1000;
xl = [-160 -130 -100]'*(pi/180);
xu = [160 130 100]'*(pi/180);

p = @(theta_1,theta_2,theta_3) [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1) + a_3*cos(theta_1 + theta_2 + theta_3); ...
                                a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1) + a_3*sin(theta_1 + theta_2 + theta_3)];


f = @(q) (1/4 * sum((pd-p(q(1),q(2),q(3)))).^2) + beta*( ~((xl(1)<q(1))&&(q(1)<xu(1))).*1 + ~((xl(2)<q(2))&&(q(2)<xu(2))).*1 + ~((xl(3)<q(3))&&(q(3)<xu(3))).*1);
%fp = @(x,xl,xu) f(x) + beta*Penalty(x,xl,xu);

G = 150;
N = 50;
D = 3;

F = 0.6;
CR = 0.9;

x = zeros(D,N);
fitness = zeros(1,N);
f_plot = zeros(1,G);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = f(x(:,i));
end 

for n=1:G
    for i=1:N
        %Mutación
        r1 = i;
        while r1 == i
            r1 = randi([1,N]);
        end

        r2 = r1;
        while r2 == r1 || r2 == i
            r2 = randi([1,N]);
        end

        r3 = r2;
        while r3 == r2 || r3 == r1 || r3 == i
            r3 = randi([1,N]);
        end

        v = x(:,r1) + F*(x(:,r2)-x(:,r3));
        %Recombinación
        u = zeros(D,1);

        for j = 1:D
            r = rand;
            if r <=  CR
                u(j) = v(j);
            else
                u(j) = x(j,i);
            end
        end

        %selección
        fitness_u = f(x(1));

        if fitness_u < fitness(i)
            x(:,i) = u;
            fitness(i) = fitness_u;
        end
    end

    f_plot(n) = min(fitness);
end

[~,igb] = min(fitness);

q = x(:,igb);
Dibujar_Manipulador(q);
plot(pd(1),pd(2),'go','LineWidth',2,'MarkerSize',10)

figure 
hold on
grid on
plot(f_plot,'b-','LineWidth',2)
title('Gráfica de convergencia')
xlabel('Iteración')
ylabel('f(x)')

function z = Penalty(x,xl,xu)
    z = 0;
    D = numel(x);

    for j=1:D
        if xl(j) < x(j)
           z = z + 0;
        else
           z = z + 1;
        end
    
        if x(j) < xu(j) 
            z = z + 0;
        else
            z = z + 1;
        end
    end
end
   