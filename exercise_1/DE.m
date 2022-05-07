clear all
close all
clc

%Sistema 1
% A = [2 3; 5 1];
% b = [7; 11];
% D = 2;
% xl = [-10; -10];
% xu = [10; 10];

%Sistema 2
% A = [2 3 4; 1 2 3; 5 1 0; 3 4 1];
% b = [19; 13; 11; 13];
% D = 3;
% xl = [-10; -10; -10];
% xu = [10; 10; 10];

%Sistema 3
% A = [2 3; 5 4; 2 5; 4 1; 0.5 0.5];
% b = [-5; 5; -15; 15; 0];
% D = 2;
% xl = [-10; -10];
% xu = [10; 10];

%Sistema 4
A = [1 -2 2 -3; 3 4 -1 1; 2 -3 2 -1; 1 1 -3 -2];
b = [15; -6; 17; -7];
D = 4;
xl = [-10; -10; -10; -10];
xu = [10; 10; 10; 10];

m = size(b,1);
f = @(x) (1/(2*m))*sum((b-A*x).^2);

G = 200;
N = 50;
F = 0.6; % o 1.2
CR = 0.9; % o 0.6

x = zeros(D,N);
fit = zeros(1,N);

for i=1:N
    x(:,i) = xl + (xu-xl).*rand(D,1);
    fit(i) = f(x(:,i));
end

f_plot = zeros(1,G);

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

    f_plot(n) = min(fit);
end

[~,igb] = min(fit);

disp(['f(x)= ' num2str(fit(igb))])
disp('x= ')
disp(x(:,igb))

%validation
disp('ans=')
disp(A*x(:,igb))
disp('b= ')
disp(b)