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

F = 0.6; % o 1.2
CR = 0.9; % o 0.6

x = zeros(D,N);
fit = zeros(1,N);

for i=1:N
    x(:,i) = xl + (xu-xl).*rand(D,1);
    fit(i) = f(x(1,i),x(2,i));
end

f_plot = zeros(1,G);

for n=1:G
    Plot_Contour(f,x,xl,xu);

    for i=1:N
        %mutación
        r1 = i;
        while r1 == i
            r1 = randi([1,N]);
        end
        
        r2 = r1;
        while r2==r1 || r2==i
            r2 = randi([1,N]);
        end
    
        r3 = r2;
        while r3==r2 || r3==r1 || r3==i
            r3 = randi([1,N]);
        end

        v = x(:,r1) + F*(x(:,r2)-x(:,r3));

        %Recombinación
        u = zeros(D,1);

        for j=1:D
            r = rand;
            if r<=CR
                u(j) = v(j);
            else
                u(j) = x(j,i);
            end
        end

        %Selección
        fit_u = f(u(1),u(2));

        if fit_u < fit(i)
            x(:,i) = u;
            fit(i) = fit_u;
        end
    end

    f_plot(n) = min(fit);
end

[~,igb] = min(fit);

Plot_Surf(f,x,xl,xu)
disp(['Mínimo global en: x=' num2str(x(1,igb)) ' y=' num2str(x(2,igb)) ', f(x,y)=' num2str(f(x(1,igb),x(2,igb)))])

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