close all
clc

f = @(x,y) x.*exp(-x.^2-y.^2);
xl = [-2, -2]';
xu = [2 2]';

% f = @(x,y) (x-2).^2+(y-2).^2;
% xl = [-5 -5]';
% xu = [5 5]';

G = 1000;    %generaciones
D = 2;      %dimension
N = 100;    %poblacion
E = 10;     %elitistas

x = zeros(D,N);
aptitud = zeros(1,N);
pm = 0.01;

%inicializacion
for i=1:N
    x(:,i) = xl + (xu-xl).*rand(D,1);
end

f_plot = zeros(1,G);

for g=1:G
    %Plot_Contour(f,x,xl,xu)
    
    %calculo aptitudes
    for i=1:N
        fx = f(x(1,i),x(2,i));
    
        if fx >= 0
            aptitud(i) = 1 / (1+fx);
        else
            aptitud(i) = 1 + abs(fx);
        end
    end

    [~,indice] = max(aptitud);
    xb = x(:,indice);
    fxb = f(xb(1),xb(2));
    f_plot(g) = fxb;

    y = zeros(D,N-E);

    for i=1:2:N-E
        r1 = Rank(aptitud); %Ruleta, Torneo
        r2 = r1;

        while r1 == r2
            r2 = Rank(aptitud); %Ruleta, Torneo
        end

        x1 = x(:,r1);
        x2 = x(:,r2);
        pc = randi([1 D]);

        alpha = rand();
        y1 = alpha*x1 + (1-alpha)*x2;
        y2 = (1-alpha)*x1 + alpha*x2;

        %y1 = [x1(1:pc) x2(pc+1:D)];
        %y2 = [x2(1:pc) x1(pc+1:D)];

        y(:,i) = y1;
        y(:,i+1) = y2;
    end

    %mutacion
    for i=1:N-E
        for j=1:D
            if rand < pm
                %y(j,i) = y(j,i) + normrnd(0,1);
                y(j,i) = xl(j) + (xu(j)-xl(j))*rand();
            end
        end
    end

    [~,I] = sort(aptitud,'descend');
    x = x(:,I);
    x = [x(:,1:E) y];
end

for i=1:N
    fx = f(x(1,i),x(2,i));

    if fx >= 0
        aptitud(i) = 1 / (1+fx);
    else
        aptitud(i) = 1 + abs(fx);
    end
end

[~,indice] = max(aptitud);
xb = x(:,indice);
fxb = f(xb(1),xb(2));

Plot_Surf(f, xb)
disp(['Mínimo global en: x=' num2str(xb(1)) ' y=' num2str(xb(2)) ', f(x,y)=' num2str(fxb)])

figure 
hold on
grid on
plot(f_plot,'b-','LineWidth',2)
title('Gráfica de convergencia')
xlabel('Iteración')
ylabel('f(x)')
return

function n = Rank (aptitud)
    [~,I] = sort(aptitud,'descend');
    N = numel(aptitud);

    rank = N:-1:1;
    rank_total = sum(rank);
    p = rank/rank_total;

    r = rand;
    p_sum = 0;

    for i=1:N
        p_sum = p_sum + p(i);

        if p_sum >= r
            n = I(i);
            return
        end
    end

    n = I(N);
end

function [n] = Torneo (aptitud)
    N = numel(aptitud);
    tao = round(N*0.3);
    I = randi(N,[1 tao]);
    [~,i] = max(aptitud(I));
    n = I(i);
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
