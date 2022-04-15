clear all
close all
clc

f = @(x,y) sin(x+y) + (x-y).^2-1.5*x+2.5*y+1; %McCormick

fp = @(x,xl,xu) f(x(1),x(2)) + 1000*Penalty(x,xl,xu);
xl = [-1.5 -3]';
xu = [4 4]';

G = 150;
N = 50;
D = 2;

F = 0.6;
CR = 0.9;

x = zeros(D,N);
fitness = zeros(1,N);
f_plot = zeros(1,G);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = f(x(1,i),x(2,i));
end 

for n=1:G
    Plot_Contour(f,x,xl,xu)

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
        fitness_u = fp(u,xl,xu);

        if fitness_u < fitness(i)
            x(:,i) = u; 
            fitness(i) = fitness_u;
        end
    end

    f_plot(n) = min(fitness);
end

[~,igb] = min(fitness);

Plot_Surf(f,x,xl,xu)
disp(['Mínimo global en: x=' num2str(x(1,igb)) ' y=' num2str(x(2,igb)) ', f(x,y)=' num2str(fitness(igb))])

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

function Plot_Contour (f,x,xl,xu)
    cla
    hold on
    grid on
    
    x_lim = linspace(xl(1),xu(1),50);
    y_lim = linspace(xl(2),xu(2),50);
    [X,Y] = meshgrid(x_lim,y_lim);
    Z = f(X,Y);

    contour(X,Y,Z,20);
    plot(x(1,:),x(2,:),'xb','LineWidth',2,'MarkerSize',10);
    plot(x(1,:),x(2,:),'or','LineWidth',2,'MarkerSize',10);

    xlabel('x','FontSize',15)
    ylabel('y','FontSize',15)
    axis([xl(1) xu(1) xl(2) xu(2)])
    pause(0.1)
end

function Plot_Surf(f,point,xl,xu)
    close all
    clc

    x_lim = linspace(xl(1),xu(1),50);
    y_lim = linspace(xl(2),xu(2),50);
    [x,y] = meshgrid(x_lim,y_lim);
    z = f(x,y);
    
    figure
    hold on
    grid on
    surf(x,y,z)
    plot3(point(1),point(2),f(point(1),point(2)),'r*','LineWidth',2,'MarkerSize',10)
    legend({'función','óptimo'},'FontSize',15)
    title('Gráfica en 3D','FontSize',15)
    xlabel('x','FontSize',15)
    ylabel('y','FontSize',15)
    zlabel('f(x,y)','FontSize',15)
    view([-20,60])
    
    figure
    hold on
    grid on
    contour(x,y,z,20)
    plot(point(1),point(2),'r*','LineWidth',2,'MarkerSize',10)
    legend({'función','óptimo'},'FontSize',15)
    title('Gráfica en 2D','FontSize',15)
    xlabel('x','FontSize',15)
    ylabel('y','FontSize',15)
end