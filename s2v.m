function vec = s2v(st)
    f = fieldnames(st);
    vec = [];
    for ctr = 1 : length(f)
        vec = [vec; st.(f{ctr})]; %#ok
    end
end
