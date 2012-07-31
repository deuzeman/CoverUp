function data = analyze(data)
    global opts;

    if strcmpi(opts.wipe, 'ON')
        close all;
    end
    
    data = prepare_data(data);

    hw = waitbar(0, 'Running calculation...');
    data = fitpions(data);
    data = fit_nucleons(data);

    boot_params = zeros(length(fieldnames(data.params)), opts.nboot);
    boot_scale = zeros(length(fieldnames(data.scale)), opts.nboot);
    
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
        boot_params(:, ctr) = s2v(samp.params);
        boot_scale(:, ctr)  = s2v(samp.scale);
        waitbar((ctr + 1) / (1 + opts.nboot), hw, 'Running bootstraps...');
    end
    data.sd_params = v2s(std(boot_params, 1, 2), fieldnames(data.params));
    data.sd_scale  = v2s(std(boot_scale, 1, 2), fieldnames(data.scale));
    
    close(hw);
    
    display_results(data);
end
