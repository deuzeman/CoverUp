function st = v2s(vec, names)
    for ctr = 1 : length(vec)
        st.(names{ctr}) = vec(ctr);
    end
end
