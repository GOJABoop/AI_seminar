clear all
close all
clc

f = @(x,y) (x-2).^2+(y-2).^2; %Sphere
%f = @(x,y) ((x.^2/4000)+(y.^2/4000)) - (cos(x).*cos(y/sqrt(2)))+1; %Griewank
%f = @(x,y) 10*2 + x.^2 + y.^2 - 10*cos(2*pi*x) - 10*cos(2*pi*y); % Rastrigin
%f = @(x,y) -((1+cos(12*sqrt(x.^2+y.^2)))./(0.5*(x.^2+y.^2)+2)); %Drop-Wave

xl = [-5 -5]';
xu = [5 5]';

G = 300;
N = 100;
D = 2;


c1 = 2.05;
c2 = 2.05;

phi =  c1 + c2;
K = (2)/(abs(2-phi-sqrt(phi^2-4*phi)));

x = zeros(D,N);
v = zeros(D,N);
fit = zeros(1,N);
xb = zeros(D,N);

f_plot = zeros(1,G);

for i =1:N
    x(:,i) = xl + (xu-xl).*rand(D,1);
    v(:,i) = randn(D,1);
    xb(:,i) = x(:,i);
    fit(i) = f(x(1,i),x(2,i));
end

for g=1:G
    Plot_Contour(f,x,xl,xu);

    for i=1:N
        fx = f(x(1,i),x(2,i));
        if fx < fit(i)
            xb(:,i) = x(:,i);
            fit(i) = fx;
        end
    end
    [f_plot(g),ig] = min(fit);
    for i=1:N
        v(:,i) = K*v(:,i) + rand() * c1 * (xb(:,i)-x(:,i)) + rand() * c2 * (xb(:,ig)-x(:,i));
        x(:,i) = x(:,i) + v(:,i);
    end
end

Plot_Surf(f, xb,xl,xu)
disp(['Mínimo global en: x=' num2str(xb(1,ig)) ' y=' num2str(xb(2,ig)) ', f(x,y)=' num2str(f(xb(1,ig),xb(2,ig)))])

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