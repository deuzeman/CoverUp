function [params, scale] = fitpions(raw_input, start, start_scale)
    global opts;

    if strcmpi(opts.wipe, 'ON')
        close all;
    end
    
    data = analyze_meta(raw_input);
    
    [params, data] =  get_starting_values(data);

    if nargin > 1
        params = start;
    end
    
    if nargin > 2
        data.scale = start_scale;
    end

    % Using these estimated scales, set up a manual iteration
    for sc_iter = 1 : 50
        % Produce a physical data set from the raw format
        data.mu  = raw_input.a_mu / (data.scale.a * data.scale.zp);               
        data.L = raw_input.L * data.scale.a;

        data.fps    = raw_input.a_fps / data.scale.a;
        data.sd_fps = raw_input.sd_a_fps / data.scale.a;
        data.mps    = raw_input.a_mps / data.scale.a;
        data.sd_mps = raw_input.sd_a_mps / data.scale.a;
        data.mps_n    = raw_input.a_mps_n / data.scale.a;
        data.sd_mps_n = raw_input.sd_a_mps_n / data.scale.a;
        
        % Perform the fit
        params = lsq(params, data);
        
        % Recalibrate the scale
        old_scale = data.scale;
        [data.scale, params] = set_scale(params, data);
    
        if sqrt(((data.scale.a - old_scale.a)^2 + (data.scale.zp - old_scale.zp)^2)) < 1e-8
            % Note that one would need to retune the scales of the data if this convergence criterion is set coarsely!
            break;
        end
    end
    
    scale = data.scale;
    display_results(params, data);
end

function [params, data] = get_starting_values(data)
    global almanac
    global opts
    
    params.f0 = 120;
    params.B0 = almanac.mpi^2 / (2 * 3.0);
    
    if strcmpi(opts.fvol, 'CDH') || strcmpi(opts.fvol, 'CWW')
        params.l1 = almanac.l1.ave;
        params.l2 = almanac.l2.ave;
    end
    
    params.l3 = almanac.l3.ave;
    params.l4 = almanac.l4.ave;    
    
    if data.meta.needs_zeta
        params.zeta = 0;
    end
    
    if data.meta.has_iso
        params.Xi3  = 0;
    end
    
    if data.meta.has_asq
        params.Dm = 0;
        params.Df = 0;
    end
    
    if data.meta.needs_Dn
        params.Dn = 0;
    end
   
    if data.meta.has_nnlo
        params.km = 0;
        params.kf = 0;
    end
    
    for idx = 1 : data.meta.num_betas - 1
        % Note that the lattice spacings are given in fm,
        % so one needs to invert them for MeV.
        params.(data.meta.fn_afac{idx}) = data.meta.a(end) / data.meta.a(idx);
        params.(data.meta.fn_zfac{idx}) = data.meta.zfac(idx);
    end
    
    data.scale.a = data.meta.a(end);
    data.scale.zp = data.meta.zp(end);
end

function data = analyze_meta(raw_input)
    global almanac
    global opts
    
    data.meta.betas = unique(raw_input.beta);
    data.meta.num_betas = length(data.meta.betas);
    
    data.meta.is_dummy = 0;
    data.meta.has_asq = strcmpi(opts.asq, 'ON');
    data.meta.has_iso = strcmpi(opts.iso, 'ON');
    data.meta.has_nnlo = strcmpi(opts.nnlo, 'ON');
    data.meta.needs_Dn = data.meta.has_asq && data.meta.has_iso;
    data.meta.needs_zeta = data.meta.has_iso || strcmpi(opts.fvol, 'CWW');
    data.meta.mps_n_mask = ~isnan(raw_input.a_mps_n);
    
    data.meta.indices = cell(size(data.meta.betas));
    data.meta.a = ones(size(data.meta.betas));
    data.meta.zp = ones(size(data.meta.betas));
    data.meta.sd_zp = ones(size(data.meta.betas));
    
    for idx = 1 : data.meta.num_betas
        data.meta.indices{idx} = find(raw_input.beta == data.meta.betas(idx));
        data.meta.a(idx) = unique(raw_input.a(data.meta.indices{idx})) / almanac.mev_fm;
        data.meta.zp(idx) = unique(raw_input.Zp(data.meta.indices{idx}));
        data.meta.sd_zp(idx) = unique(raw_input.sd_Zp(data.meta.indices{idx}));
    end
    
    if data.meta.num_betas < 2
        return
    end
    
    data.meta.fn_zfac = cell(data.meta.num_betas - 1, 1);
    data.meta.fn_afac = cell(data.meta.num_betas - 1, 1);
    data.meta.zfac = ones(data.meta.num_betas - 1, 1);
    data.meta.sd_zfac = zeros(data.meta.num_betas - 1, 1);
    
    for idx = 1 : data.meta.num_betas - 1
        data.meta.fn_zfac{idx} = sprintf('z_%u_over_%u', 1000 * data.meta.betas(idx), 1000 * data.meta.betas(end));
        data.meta.fn_afac{idx} = sprintf('a_%u_over_%u', 1000 * data.meta.betas(idx), 1000 * data.meta.betas(end));
        data.meta.zfac(idx) = data.meta.zp(idx) / data.meta.zp(end);
        data.meta.sd_zfac(idx) = data.meta.zfac(idx) * ...
            sqrt((data.meta.sd_zp(idx) / data.meta.zp(idx))^2 + (data.meta.sd_zp(end) / data.meta.zp(end))^2);
    end
end
