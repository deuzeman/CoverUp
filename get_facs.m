function data = get_facs(data)
    data.afac = ones(size(data.mu));
    data.zfac = ones(size(data.mu));
   
    for idx = 1 : data.meta.num_betas - 1
        data.afac(data.meta.indices{idx}) = data.params.(data.meta.fn_afac{idx});
        data.zfac(data.meta.indices{idx}) = data.params.(data.meta.fn_zfac{idx});
    end
end
