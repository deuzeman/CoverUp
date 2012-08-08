function [result, data] = chi(params, data)
    global opts;
   
    data.params = params;
    data = calculate_predictions(data);

    % Calculate the deviation in the pion mass and decay constant
    data.dev.mps   = data.mps   - sqrt((data.inf_mps2   + data.asq_mps2)   .* data.fvol_mps2);
    data.dev.fps   = data.fps   -      (data.inf_fps    + data.asq_fps)    .* data.fvol_fps;
    data.dev.mps_n = data.mps_n - sqrt((data.inf_mps2_n + data.asq_mps2_n) .* data.fvol_mps2_n);
    
    data.chi.mps = data.dev.mps ./ data.sd_mps;
    data.chi.fps = data.dev.fps ./ data.sd_fps;
    if data.meta.has_iso
        data.chi.mps_n = data.dev.mps_n ./ data.sd_mps_n;
        data.chi.mps_n(isnan(data.chi.mps_n)) = [];
    end
    % Calculate the deviations from the priors, if requested   
    if ~strcmpi(opts.priors, 'OFF')
        data = priors_chisq(data);
    end

    result = data.chi;
end

function data = priors_chisq(data)
    global almanac;
    global opts;
       
    has_l12 = strcmpi(opts.fvol, 'CDH') || strcmpi(opts.fvol, 'CWW');
    
    num_priors = 2 * (data.meta.num_betas - 1) + 2 * has_l12 + 2 * strcmpi(opts.priors, 'ON');

    data.chi.priors = zeros(num_priors, 1);
    
    idx = 1;
    for beta_ctr = 1 : data.meta.num_betas - 1
        if ~strcmpi(opts.priors, 'NOZP')
            data.chi.priors(idx) = (data.params.(data.meta.fn_zfac{ceil(idx / 2)}) - data.meta.zfac(ceil(idx / 2))) / data.meta.sd_zfac(ceil(idx / 2)); 
        end
        data.chi.priors(idx + 1) = (data.params.(data.meta.fn_afac{ceil(idx / 2)}) - data.meta.afac(ceil(idx / 2))) / data.meta.sd_afac(ceil(idx / 2)); 
        idx = idx + 2;
    end
    
    if has_l12
        data.chi.priors(idx)     = (data.params.l1 - almanac.l1.ave) / almanac.l1.std;
        data.chi.priors(idx + 1) = (data.params.l2 - almanac.l2.ave) / almanac.l2.std;
        idx = idx + 2;
    end
    
    if strcmpi(opts.priors, 'ON')
        data.chi.priors(idx)     = (data.params.l3 - almanac.l3.ave) / almanac.l3.std;
        data.chi.priors(idx + 1) = (data.params.l4 - almanac.l4.ave) / almanac.l4.std;
    end
end
