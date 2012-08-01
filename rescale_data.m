function data = rescale_data(data)
    if data.meta.is_dummy
        data.afac = ones(size(data.mu));
        data.zfac = ones(size(data.mu));
        return
    end

    data = get_facs(data);
    data.scale.a2 = data.afac.^2 .* data.scale.a.^2;

    data.mps = data.mps ./ data.afac;
    data.sd_mps = data.sd_mps ./ data.afac;
    data.fps = data.fps ./ data.afac;
    data.sd_fps = data.sd_fps ./ data.afac;
    data.mu = data.mu ./ (data.afac .* data.zfac);
    data.L = data.L .* data.afac;

    data.mps_n = data.mps_n ./ data.afac;
    data.sd_mps_n = data.sd_mps_n ./ data.afac;
    
    data.mk = data.mk ./ data.afac;
    data.sd_mk = data.sd_mk ./ data.afac;
    
    data.md = data.md ./ data.afac;
    data.sd_md = data.sd_md ./ data.afac;
    
    data.mn = data.mn ./ data.afac;
    data.sd_mn = data.sd_mn ./ data.afac;
end
