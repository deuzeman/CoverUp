function data = lsq(data)
    names = fieldnames(data.params);
    problem.solver = 'lsqnonlin';
    problem.options = optimset('MaxIter', 200, 'Display', 'off');
    problem.objective = @(x)(s2v(chi(v2s(x, names), data)));
    problem.x0 = s2v(data.params);
    vres = lsqnonlin(problem);
    data.params = v2s(vres, names);
end
