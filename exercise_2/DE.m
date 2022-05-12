clear all
close all
clc

data = readtable('data_lineal_1.csv'); % 1, 2, 3

X = data.x;
Y = data.y;

n = size(X,1);
X = [ones(n,1) X];

f = @(w) (1/(2*n))*sum((Y-(X*w)).^2);

xl = [-10; -10];
xu = [10; 10];

G = 150;
N = 50;
D = 2;

F = 0.6; % o 1.2
CR = 0.9; % o 0.6

x = zeros(D,N);
fit = zeros(1,N);

for i=1:N
    x(:,i) = xl + (xu-xl).*rand(D,1);
    fit(i) = f(x(:,i));
end

for n=1:G
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
        fit_u = f(u);

        if fit_u < fit(i)
            x(:,i) = u;
            fit(i) = fit_u;
        end
    end
end

[~,igb] = min(fit);

w = x(:,igb);
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
