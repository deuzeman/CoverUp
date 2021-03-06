function data = prepare_data(raw)
    global almanac
    global opts
    
    data.raw = raw;
    
    data.meta.betas = unique(raw.beta);
    data.meta.num_betas = length(data.meta.betas);

    data.meta.is_dummy = 0;
    data.meta.has_asq = strcmpi(opts.asq, 'ON');
    data.meta.has_iso = strcmpi(opts.iso, 'ON');
    data.meta.has_nnlo = strcmpi(opts.nnlo, 'ON');
    data.meta.use_corr = strcmpi(opts.corr, 'ON');
    data.meta.needs_Dn = data.meta.has_asq && data.meta.has_iso;
    data.meta.needs_zeta = data.meta.has_iso || strcmpi(opts.fvol, 'CWW');
    data.meta.has_l12 = strcmpi(opts.fvol, 'CDH') || strcmpi(opts.fvol, 'CWW');
    data.meta.mps_n_mask = ~isnan(raw.a_mps_n);
    
    data.meta.indices = cell(size(data.meta.betas));
    data.meta.a = ones(size(data.meta.betas));
    data.meta.zp = ones(size(data.meta.betas));
    data.meta.sd_zp = ones(size(data.meta.betas));
    
    for idx = 1 : data.meta.num_betas
        data.meta.indices{idx} = find(raw.beta == data.meta.betas(idx));
        data.meta.a(idx) = unique(raw.a(data.meta.indices{idx})) / almanac.mev_fm;
        data.meta.zp(idx) = unique(raw.Zp(data.meta.indices{idx}));
        data.meta.sd_zp(idx) = unique(raw.sd_Zp(data.meta.indices{idx}));
        data.meta.r0(idx) = unique(raw.r0_a(data.meta.indices{idx}));
        data.meta.sd_r0(idx) = unique(raw.sd_r0_a(data.meta.indices{idx}));
    end
    
    data.scale.a = data.meta.a(end);
    data.scale.zp = data.meta.zp(end);
   
    % Plug in some reasonable starting values
    data.params.f0 = 120;
    data.params.B0 = almanac.mpi^2 / (2 * 3.0);
    
    if strcmpi(opts.fvol, 'CDH') || strcmpi(opts.fvol, 'CWW')
        data.params.l1 = almanac.l1.ave;
        data.params.l2 = almanac.l2.ave;
    end
    
    data.params.l3 = almanac.l3.ave;
    data.params.l4 = almanac.l4.ave;    
    
    if data.meta.needs_zeta
        data.params.zeta = 0;
    end
    
    if data.meta.has_iso
        data.params.Xi3  = 0;
    end
    
    if data.meta.has_asq
        data.params.Dm = 0;
        data.params.Df = 0;
    end
    
    if data.meta.needs_Dn
        data.params.Dn = 0;
    end
   
    if data.meta.has_nnlo
        data.params.km = 0;
        data.params.kf = 0;
    end
    
    if data.meta.num_betas < 2
        return
    end
    
    data.meta.fn_a = cell(data.meta.num_betas, 1);
    data.meta.fn_zfac = cell(data.meta.num_betas - 1, 1);
    data.meta.fn_afac = cell(data.meta.num_betas - 1, 1);
    data.meta.zfac = ones(data.meta.num_betas - 1, 1);
    data.meta.sd_zfac = zeros(data.meta.num_betas - 1, 1);
    data.meta.afac = ones(data.meta.num_betas - 1, 1);
    data.meta.sd_afac = zeros(data.meta.num_betas - 1, 1);

    for idx = 1 : data.meta.num_betas
        data.meta.fn_a{idx} = sprintf('a_%u', 1000 * data.meta.betas(idx));
    end
    
    for idx = 1 : data.meta.num_betas - 1
        data.meta.fn_zfac{idx} = sprintf('z_%u_over_%u', 1000 * data.meta.betas(end), 1000 * data.meta.betas(idx));
        data.meta.fn_afac{idx} = sprintf('a_%u_over_%u', 1000 * data.meta.betas(end), 1000 * data.meta.betas(idx));
        data.meta.zfac(idx) = data.meta.zp(idx) / data.meta.zp(end);
        data.meta.sd_zfac(idx) = data.meta.zfac(idx) * ...
            sqrt((data.meta.sd_zp(idx) / data.meta.zp(idx))^2 + (data.meta.sd_zp(end) / data.meta.zp(end))^2);
        data.meta.afac(idx) = data.meta.r0(end) / data.meta.r0(idx);
        data.meta.sd_afac(idx) = data.meta.afac(idx) * ...
            sqrt((data.meta.sd_r0(idx) / data.meta.r0(idx))^2 + (data.meta.sd_r0(end) / data.meta.r0(end))^2);
    end
    
    for idx = 1 : data.meta.num_betas - 1
        % Note that the lattice spacings are given in fm,
        % so one needs to invert them for MeV.
        data.params.(data.meta.fn_afac{idx}) = data.meta.r0(end) / data.meta.r0(idx);
        data.params.(data.meta.fn_zfac{idx}) = data.meta.zfac(idx);
    end
end
