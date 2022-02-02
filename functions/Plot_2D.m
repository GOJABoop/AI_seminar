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