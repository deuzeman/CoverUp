function result = lsq(params, data)
    names = fieldnames(params);
    problem.solver = 'lsqnonlin';
    problem.options = optimset('MaxIter', 200, 'Display', 'off');
    problem.objective = @(x)(s2v(chi(v2s(x, names), data)));
    problem.x0 = s2v(params);
    vres = lsqnonlin(problem);
    result = v2s(vres, names);
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
