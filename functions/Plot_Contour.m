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