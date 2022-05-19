clear all
close all
clc

a_1 = 0.35;
a_2 = 0.35;
a_3 = 0.25;

beta=1000;
xl=[-160 -130 -100]'*(pi/180);
xu=[160 130 100]'*(pi/180);

pd = [-0.5; 0.5];

p = @(theta_1,theta_2,theta_3) [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1) + a_3*cos(theta_1 + theta_2 + theta_3); ...
                                a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1) + a_3*sin(theta_1 + theta_2 + theta_3)];

f = @(q) (1/4)*sum((pd-p(q(1),q(2),q(3))).^2) + beta*(~((xl(1)<q(1))&&(q(1)<xu(1))).*1 + ~((xl(2)<q(2))&&(q(2)<xu(2))).*1 + ~((xl(3)<q(3))&&(q(3)<xu(3))).*1);

q = TLBO(f,xu,xl);
disp(['q(1)=' num2str(q(1)) ' q(2)='  num2str(q(2)) ', q(3)='  num2str(q(3))])

figure
cla
Dibujar_Manipulador(q)
plot(pd(1),pd(2),'go','LineWidth',2,'MarkerSize',10)


function Dibujar_Manipulador (q)
    a_1 = 0.35;
    a_2 = 0.35;
    a_3 = 0.25;
    
    theta_1 = q(1);
    theta_2 = q(2);
    theta_3 = q(3);
    
    t01 = [a_1*cos(theta_1); a_1*sin(theta_1)];
	
    t02 = [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1); ...
           a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1)];
    
    t03 = [a_2*cos(theta_1 + theta_2) + a_1*cos(theta_1) + a_3*cos(theta_1 + theta_2 + theta_3); ...
           a_2*sin(theta_1 + theta_2) + a_1*sin(theta_1) + a_3*sin(theta_1 + theta_2 + theta_3)];
   
    hold on 
    grid on

    line([0 t01(1)],[0 t01(2)],'color',[0 0 0],'LineWidth',4)
    line([t01(1) t02(1)],[t01(2) t02(2)],'color',[0 0 0],'LineWidth',4)
    line([t02(1) t03(1)],[t02(2) t03(2)],'color',[0 0 0],'LineWidth',4)

    plot(0,0,'ko','MarkerSize',15,'LineWidth',2)
    plot(t01(1),t01(2),'ko','MarkerSize',15,'LineWidth',2)
    plot(t02(1),t02(2),'ko','MarkerSize',15,'LineWidth',2)
    
    axis([-1 1 -1 1])
end

function Plot_F(f_plot)
    figure
    cla
    hold on
    grid on
    plot(f_plot, 'b-', 'LineWidth', 2)
    title('Grafica de convergencia')
    xlabel('Iteracion')
    ylabel('f(x)')
end

function q = TLBO(f,xu,xl) %Optimización Basada en la Enseñanza-Aprendizaje
    D = 3;
    N = 50;
    G = 100;
    
    x = zeros(D,N);
    fitness = zeros(1,N);
    f_plot=zeros(1,G);
    
    for i=1:N
        x(:,i) = xl+(xu-xl).*rand(D,1);
        fitness(i) = f(x);
    end
    
    for g=1:G
        
        for i=1:N
            [~,t] = min(fitness);
            Tf = randi([1 2]);
            c = zeros(D,1);

            for j=1:D
                X = mean(x(j,:));
                c(j) = x(j,i) + rand()*(x(j,t)-Tf*X);
            end
            
            fc = f(c);
            
            if fc<fitness(i)
                x(:,i) = c;
                fitness(i) = fc;
            end
            
            c = zeros(D,1);
            k = i;
            while k==i
                k = randi([1 N]);
            end
            
            if fitness(i)<fitness(k)
                for j=1:D
                    c(j) = x(j,i) + rand()*(x(j,i)-x(j,k));
                end
            else
                for j=1:D
                    c(j) = x(j,i) + rand()*(x(j,k)-x(j,i));
                end
            end
            
            fc = f(c);
            
            if fc<fitness(i)
                x(:,i) = c;
                fitness(i) = fc;
            end
        end
        [~,igb] = min(fitness);
        f_plot(g)=fitness(igb);
    end 

    q=x;
    Plot_F(f_plot);
end