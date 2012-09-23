function [result, data] = chi(params, data)
    global opts;
   
    data.params = params;
    data = calculate_predictions(data);

    % Calculate the deviation in the pion mass and decay constant
    data.dev.mps =   data.mps - sqrt((data.inf_mps2 + data.asq_mps2) .* data.fvol_mps2);
    data.dev.fps =   data.fps - (data.inf_fps  + data.asq_fps)  .* data.fvol_fps;   
    data.dev.mps_n = data.mps_n - sqrt((data.inf_mps2_n + data.asq_mps2_n) .* data.fvol_mps2_n);
    
    cov_mat = zeros(length(data.sd_mps) + length(data.sd_fps));
    for ctr = 1 : length(data.sd_mps)
        cov_mat(2 * ctr - 1, 2 * ctr - 1) = data.sd_mps(ctr)^2;
        cov_mat(2 * ctr    , 2 * ctr    ) = data.sd_fps(ctr)^2;
        if data.meta.use_corr
            cov_mat(2 * ctr - 1, 2 * ctr    ) = data.raw.corr(ctr) * data.sd_mps(ctr) * data.sd_fps(ctr);
            cov_mat(2 * ctr    , 2 * ctr - 1) = data.raw.corr(ctr) * data.sd_fps(ctr) * data.sd_mps(ctr);
        end
    end
    [cma, cmb, cmc] = svd(cov_mat);
    weight = cma * diag(1 ./ sqrt(diag(cmb))) * cmc';
        
    chi_tot = [data.dev.mps, data.dev.fps]';
    chi_tot = reshape(chi_tot, numel(chi_tot), 1);
    chi_tot = weight * chi_tot;
    chi_tot = reshape(chi_tot, 2, numel(chi_tot) / 2)';
    
    data.chi.mps = chi_tot(:, 1);
    data.chi.fps = chi_tot(:, 2);

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
        
    num_priors = (data.meta.num_betas - 1) + 2 * data.meta.has_l12 + (data.meta.num_betas - 1) * ~(strcmpi(opts.priors, 'NOZP') || strcmpi(opts.priors, 'NOR0')) + 2 * strcmpi(opts.priors, 'ON');

    data.chi.priors = zeros(num_priors, 1);
    
    idx = 1;
    for beta_ctr = 1 : data.meta.num_betas - 1
        if ~strcmpi(opts.priors, 'NOZP')
            data.chi.priors(idx) = (data.params.(data.meta.fn_zfac{beta_ctr}) - data.meta.zfac(beta_ctr)) / data.meta.sd_zfac(beta_ctr); 
            idx = idx + 1;
        end
        if ~strcmpi(opts.priors, 'NOR0')
            data.chi.priors(idx) = (data.params.(data.meta.fn_afac{beta_ctr}) - data.meta.afac(beta_ctr)) / data.meta.sd_afac(beta_ctr); 
            idx = idx + 1;
        end
    end
    
    if data.meta.has_l12
        data.chi.priors(idx)     = (data.params.l1 - almanac.l1.ave) / almanac.l1.std;
        data.chi.priors(idx + 1) = (data.params.l2 - almanac.l2.ave) / almanac.l2.std;
        idx = idx + 2;
    end
    
    if strcmpi(opts.priors, 'ON')
        data.chi.priors(idx)     = (data.params.l3 - almanac.l3.ave) / almanac.l3.std;
        data.chi.priors(idx + 1) = (data.params.l4 - almanac.l4.ave) / almanac.l4.std;
    end
end
