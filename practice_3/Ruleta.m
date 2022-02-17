function indice = Ruleta(aptitup)
    p = aptitup/sum(aptitup);
    N = numel(aptitup);
    r = rand();
    p_sum = 0;

    for i=1:N
        p_sum = p_sum + p(i);

        if p_sum >= r
            indice = i;
            return
        end
    end

    indice = N;
end

