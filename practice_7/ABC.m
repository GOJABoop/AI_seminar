clear all
close all
clc

f = @(x,y) (x-2).^2+(y-2).^2; %Sphere
%f = @(x,y) ((x.^2/4000)+(y.^2/4000)) - (cos(x).*cos(y/sqrt(2)))+1; %Griewank
%f = @(x,y) 10*2 + x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y); % Rastrigin
%f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2)); %Drop-Wave

xl = [-5 -5]';
xu = [5 5]';

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
    fitness(i) = f(x(1,i),x(2,i));
end

for g=1:G
    %Plot_Contour(f,x,xl,xu)
    
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
        fv = f(v(1),v(2));

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
        fv = f(v(1),v(2));

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
            fitness(i) = f(x(1,i),x(2,i));
            l(i) = 0;
        end
    end

    f_plot(g) = min(fitness);
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