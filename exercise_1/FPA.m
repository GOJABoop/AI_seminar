clear all
close all
clc

%Sistema 1
A = [2 3; 5 1];
b = [7; 11];
D = 2;
xl = [-10; -10];
xu = [10; 10];

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
% A = [1 -2 2 -3; 3 4 -1 1; 2 -3 2 -1; 1 1 -3 -2];
% b = [15; -6; 17; -7];
% D = 4;
% xl = [-10; -10; -10; -10];
% xu = [10; 10; 10; 10];

m = size(b,1);
f = @(x) (1/(2*m))*sum((b-A*x).^2);

G = 200;
N = 50;
p = 0.8;
lambda = 1.5;
sigma2 = (((gamma(1+lambda))/(lambda*gamma((1+lambda)/2)))*((sin((pi*lambda)/2))/(2^((lambda-1)/2))))^(1/lambda);
x = zeros(D,N);
fitness = zeros(1,N);
f_plot = zeros(1,G);

for i=1:N
    x(:,i) = xl+(xu-xl).*rand(D,1);
    fitness(i) = f(x(:,i));
end 

for t=1:G  
    [~,igb] = min(fitness);

    for i=1:N
        if rand() < p
            u = normrnd(0,sigma2,[D 1]);
            v = normrnd(0,1,[D 1]);
            L = u./(abs(v).^(1/lambda));
            y = x(:,i) + L.*(x(:,igb)-x(:,i));
        else
            j = i;
            while j==i
                j = randi([1 N]);
            end

            k = j;
            while k==j && k~=i
                k = randi([1 N]);
            end
            
            y = x(:,i) + rand()*(x(:,j)-x(:,k));
        end

        fy = f(y);

        if fy < fitness(i)
            x(:,i) = y;
            fitness(i) = fy;
        end
    end

    f_plot(t) = fitness(i);
end

[~,igb] = min(fitness);

disp(['f(x)= ' num2str(fitness(igb))])
disp('x= ')
disp(x(:,igb))

%validation
disp('ans=')
disp(A*x(:,igb))
disp('b= ')
disp(b)