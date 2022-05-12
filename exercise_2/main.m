clear all
close all
clc

data = readtable('data_lineal_1.csv'); % 1, 2, 3

X = data.x;
Y = data.y;

n = size(X,1);
X = [ones(n,1) X];

f = @(w) (1/(2*n))*sum((Y-(X*w)).^2);


w = [0.5 0.5]';
disp([' Coeficientes: w = [' num2str(w') '], f(w) = ' num2str(f(w))])


Yp = X*w;
R2 = 1 - sum((Y-Yp).^2)/sum((Y-mean(Y)).^2);
disp([' R2 score = ' num2str(R2)])


figure
hold on
grid on
plot(X(:,2),Y,'bo','LineWidth',2,'MarkerSize',10)
plot(X(:,2),Yp,'g-','LineWidth',3,'MarkerSize',10)
xlabel('x','FontSize',13)
ylabel('y','FontSize',13)
legend({'Muestras','Regresion'},'FontSize',13)
