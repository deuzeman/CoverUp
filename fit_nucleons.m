function data = fit_nucleons(data)
    global almanac;
    
    data = get_facs(data);
    names = {'mn0', 'c0', 'ga', 'Dnuc'};
    mask = ~isnan(data.mn); 
    
    A = [data.mps(mask), data.mps(mask).^2, - (3 * data.mps(mask).^3) ./ (16 * pi * almanac.fpi), data.afac(mask).^2];
    x = lscov(A, data.mn(mask), data.sd_mn(mask).^-2);
    data.params_nuc = v2s(x, names);
    data.chi_nuc = (data.mn(mask) - A * x) ./ data.sd_mn(mask);
end
