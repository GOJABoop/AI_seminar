close all
clc

%-------------------Función 1----------------------------
f = @(x,y) x.*exp(-x.^2-y.^2);
g = @(x,y) [(1-2*x^2)*exp(-x^2-y^2) -2*x*y*exp(-x^2-y^2)]';
%lower_limit = [-3 -3]';
%upper_limit = [3 3]';
xi = [-1 1]';

h = 0.1;

for i = 1:100
    %Plot_Contour(f,xi,lower_limit,upper_limit)
    xi = xi - h * g(xi(1),xi(2));
end

Plot_Surf(f,xi)

disp('Función 1')
disp(['Mínimo global en: x=' num2str(xi(1)) ', y=' num2str(xi(2)) ...
      ', f(x,y)=' num2str(f(xi(1),xi(2)))])

%-------------------Función 2----------------------------
f =  @(x,y) 2*x.^2-8*x+8;
g = @(x,y) [4*x-8 0]';
%lower_limit = [-3 -3]';
%upper_limit = [2.5 3]';
xi = [-3 0]';

h = 0.1;

for i = 1:100
    %Plot_Contour(f,xi,lower_limit,upper_limit)
    xi = xi - h * g(xi(1),xi(2));
end

Plot_Surf(f,xi)
disp('Función 2')
disp(['Mínimo global en: x=' num2str(xi(1)) ', y=' num2str(xi(2)) ...
      ', f(x,y)=' num2str(f(xi(1),xi(2)))])


function Plot_Contour (f,x,lower_limit,upper_limit)
    cla
    hold on
    grid on
    
    x_lim = linspace(lower_limit(1),upper_limit(1),50);
    y_lim = linspace(lower_limit(2),upper_limit(2),50);
    [X,Y] = meshgrid(x_lim,y_lim);
    Z = f(X,Y);

    contour(X,Y,Z,20);
    plot(x(1,:),x(2,:),'xb','LineWidth',2,'MarkerSize',10);
    plot(x(1,:),x(2,:),'or','LineWidth',2,'MarkerSize',10);

    xlabel('x','FontSize',15)
    ylabel('y','FontSize',15)

    axis([lower_limit(1) upper_limit(1) lower_limit(2) upper_limit(2)])
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

