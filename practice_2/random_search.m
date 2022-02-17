close all
clc

%f = @(x,y) x.*exp(-x.^2-y.^2);
f = @(x,y) (x-2).^2+(y-2).^2;
%f =  @(x,y) x.^2-4*x+8*y.^2-4*y;
xl = [-2, -2]';
xu = [2 2]';

G = 100000;
D = 2;

xb = [0 0]';
fb = inf;

f_plot =  zeros(1,G);

for g=1:G
    x = xl + (xu-xl).*rand(D,1);
    fx = f(x(1),x(2));
    if fx < fb
        fb = fx;
        xb = x;
    end
    %Plot_Contour(f,[x xb],xl,xu);
    f_plot(g) = fb;
end

Plot_Surf(f, xb)
disp(['Mínimo global en: x=' num2str(xb(1)) ' y=' num2str(xb(2)) ', f(x,y)=' num2str(fb)])

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

function Plot_Surf(f,point)
    close all
    clc
    
    x_lim = linspace(-5,5,50);
    y_lim = linspace(-5,5,50);
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
