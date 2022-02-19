%(mu+1)-ES
clear all
close all
clc

f = @(x,y) x.*exp(-x.^2-y.^2);
xl = [-2, -2]';
xu = [2 2]';

% f = @(x,y) (x-2).^2+(y-2).^2;
% xl = [-5 -5]';
% xu = [5 5]';

G = 500;    %generaciones
D = 2;      %dimensiones
mu = 20;      

x = zeros(D,mu+1);
fit = zeros(1,mu+1);
sigma = zeros(D,mu+1);

for i=1:mu
   x(:,i) = xl+(xu-xl).*rand(D,1);
   fit(i) = f(x(1,i),x(2,i));
   sigma(:,i) = 0.1*rand(D,1);
end

f_plot = zeros(1,G);

for t=1:G
    r1 = randi([1 mu]);
    r2 = r1;

    while r2 == r1
        r2 = randi([1 mu]);
    end
    
    x(:,mu+1) = (x(:,r1)+x(:,r2)) /2;
    sigma(:,mu+1) = (sigma(:,r1)+sigma(:,r2)) /2;
    r = normrnd(0,sigma(:,mu+1));
    x(:,mu+1) =  x(:,mu+1) + r;
    
    for i=1:mu+1
       fit(i) = f(x(1,i),x(2,i));
    end
    
    [~,I] = sort(fit);
    
    x = x(:,I);
    sigma = sigma(:,I);
    fit = fit(I);

    f_plot(t) = fit(1);
    %Plot_Contour(f,x(:,1:mu),xl,xu);
end

xb = x(:,1);
Plot_Surf(f, xb,xl,xu)
disp(['Mínimo global en: x=' num2str(xb(1)) ' y=' num2str(xb(2)) ', f(x,y)=' num2str(f(xb(1),xb(2)))])

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

