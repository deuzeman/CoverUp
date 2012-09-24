function data = analyze(data)
    global opts;

    if strcmpi(opts.wipe, 'ON')
        close all;
    end
    
    data = prepare_data(data);

    hw = waitbar(0, 'Running calculation...');
    data = fitpions(data);
    if (sum(~isnan(data.mn)) > 5)
        data = fit_nucleons(data);
    end

    data = etmc_convention(data);
    
    boot_params = zeros(length(fieldnames(data.params)), opts.nboot);
    boot_scale = zeros(length(fieldnames(data.scale)), opts.nboot);
    boot_params_etmc = zeros(length(fieldnames(data.params_etmc)), opts.nboot);
    boot_scale_etmc = zeros(length(fieldnames(data.scale_etmc)), opts.nboot);
    
    waitbar(1 / (1 + opts.nboot), hw, 'Running bootstraps...');
    for ctr = 1 : opts.nboot
        noise_fps = randn(size(data.raw.a_fps));
        noise_mps = data.raw.corr.^2 .* noise_fps + (1 - data.raw.corr).^2 .* randn(size(data.raw.a_mps));
        noise_mps_n = randn(size(data.raw.a_fps));
        samp = data;
        samp.raw.a_fps = samp.raw.a_fps + data.raw.sd_a_fps .* noise_fps;
        samp.raw.a_mps = samp.raw.a_mps + data.raw.sd_a_mps .* noise_mps;
        samp.raw.a_mps_n = samp.raw.a_mps_n + data.raw.sd_a_mps .* noise_mps_n;
        samp = fitpions(samp);
        samp = etmc_convention(samp);
        boot_params(:, ctr) = s2v(samp.params);
        boot_scale(:, ctr)  = s2v(samp.scale);
        boot_params_etmc(:, ctr) = s2v(samp.params_etmc);
        boot_scale_etmc(:, ctr) = s2v(samp.scale_etmc);
        waitbar((ctr + 1) / (1 + opts.nboot), hw, 'Running bootstraps...');
    end
      
    data.sd_params = v2s(std(boot_params, 1, 2), fieldnames(data.params));
    data.sd_scale  = v2s(std(boot_scale, 1, 2), fieldnames(data.scale));

    data.sd_params_etmc = v2s(std(boot_params_etmc, 1, 2), fieldnames(data.params_etmc));
    data.sd_scale_etmc  = v2s(std(boot_scale_etmc, 1, 2), fieldnames(data.scale_etmc));

    close(hw);
    
    data = display_results(data);
end

function data = etmc_convention(data)
    global almanac;
    
    % Convert LECs to length scales
    if data.meta.has_l12
        data.params_etmc.l1 = sqrt(exp(data.params.l1)) * almanac.mpi;
        data.params_etmc.l2 = sqrt(exp(data.params.l2)) * almanac.mpi;
    end
    data.params_etmc.l3 = sqrt(exp(data.params.l3)) * almanac.mpi;
    data.params_etmc.l4 = sqrt(exp(data.params.l4)) * almanac.mpi;
    if data.meta.has_iso
        data.params_etmc.Xi3 = sqrt(exp(data.params.Xi3)) * almanac.mpi;
    end
    
    % Express all observables in terms of f0
    data.params_etmc.B0 = data.params.B0 / data.params.f0;
    if data.meta.has_l12
        data.params_etmc.l1 = data.params_etmc.l1 / data.params.f0;
        data.params_etmc.l2 = data.params_etmc.l2 / data.params.f0;
    end
    data.params_etmc.l3 = data.params_etmc.l3 / data.params.f0;
    data.params_etmc.l4 = data.params_etmc.l4 / data.params.f0;
    
   if data.meta.needs_zeta
        data.params_etmc.zeta = data.params.zeta / data.params.f0^4;
    end
    
    if data.meta.has_iso
        data.params_etmc.Xi3  = data.params_etmc.Xi3 / data.params.f0;
    end
    
    if data.meta.has_asq
        data.params_etmc.Dm = data.params.Dm / data.params.f0^2;
        data.params_etmc.Df = data.params.Df / data.params.f0^2;
    end
    
    if data.meta.needs_Dn
        data.params_etmc.Dn = data.params.Dn / data.params.f0^2;
    end
   
    if data.meta.has_nnlo
        data.params_etmc.km = data.params.km / data.params.f0^2;
        data.params_etmc.kf = data.params.kf / data.params.f0^2;
    end

    for idx = 1 : data.meta.num_betas
        data.scale_etmc.(data.meta.fn_a{idx}) = data.scale.(data.meta.fn_a{idx}) * data.params.f0;
    end
end