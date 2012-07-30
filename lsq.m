function data = lsq(data)
    names = fieldnames(data.params);
    problem.solver = 'lsqnonlin';
    problem.options = optimset('MaxIter', 200, 'Display', 'off');
    problem.objective = @(x)(s2v(chi(v2s(x, names), data)));
    problem.x0 = s2v(data.params);
    vres = lsqnonlin(problem);
    data.params = v2s(vres, names);
end

function st = v2s(vec, names)
    for ctr = 1 : length(vec)
        st.(names{ctr}) = vec(ctr);
    end
end

function vec = s2v(st)
    f = fieldnames(st);
    vec = [];
    for ctr = 1 : length(f)
        vec = [vec; st.(f{ctr})]; %#ok
    end
end
