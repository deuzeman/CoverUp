function data = fit_nucleons(data)
    global almanac;
    
    data = get_facs(data);
    names = {'mn0', 'c0', 'ga', 'Dnuc'};
    mask = ~isnan(data.mn); 
    
    A = [ones(size(data.mps(mask))), data.mps(mask).^2,  -(3 * data.mps(mask).^3) ./ (16 * pi * almanac.fpi^2), data.afac(mask).^2 .* data.scale.a.^2];
    x = lscov(A, data.mn(mask), data.sd_mn(mask).^-2);
    
    data.chi_nuc = (data.mn(mask) - A * x) ./ data.sd_mn(mask);
    
    data.params_nuc = v2s(x, names);
    data.params_nuc.c0 = 0.25 * data.params_nuc.c0;
    data.params_nuc.ga = sqrt(data.params_nuc.ga);
end
