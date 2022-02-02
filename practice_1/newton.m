close all
clc

f = @(x) 4*x.^3 - 80*x.^2 + 400*x;
fp = @(x) 12*x.^2 - 160*x + 400;
fpp = @(x) 24*x -160;

limiteInferior = 0;
limiteSuperior = 10;
xi = 5;
valorInicial = xi;
n = 10;

for i = 1 : n
   Plot_2D(fp,xi,limiteInferior,limiteSuperior) %Gráfica
   if fpp(xi) ~= 0
    xi = xi - (fp(xi)/fpp(xi));
   end
end

%Desplegar gráfica
figure
hold on
grid on

xp = limiteInferior:0.01:limiteSuperior;
plot(xp,f(xp),'b-','LineWidth',3,'MarkerSize',12)
plot(xp,fp(xp),'r-','LineWidth',2,'MarkerSize',12);
plot(xp,fpp(xp),'g:','LineWidth',2,'MarkerSize',12);
plot(xi,f(xi),'*r','LineWidth',3,'MarkerSize',12);
legend({'f(x)','f''(x)','f''''(x)','óptimo'},'FontSize',15)
title('Gráfica en 2D','FontSize',15)
xlabel('x')
ylabel('f(x)')

%Resultados
disp(['Valor inicial: ' num2str(valorInicial)])
disp(['Raíz en la f''(x): x = ' num2str(xi)])
if fpp(xi) >= 0
    disp(['Minímo en ' num2str(xi)])
else
    disp(['Máximo en ' num2str(xi)])
end


function Plot_2D(f,x, limiteInferior, limiteSuperior)
    cla
    hold on
    grid on

    xp = limiteInferior:0.01:limiteSuperior;
    plot(xp,f(xp),'b-','LineWidth',3,'MarkerSize',12)
    plot(x,f(x),'*r','LineWidth',3,'MarkerSize',12);
    xlabel('x')
    ylabel('f(x)')

    pause(0.5)
end