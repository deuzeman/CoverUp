function [scale, params] = set_scale(params, data)  
  global almanac;
  
  phys_rat = almanac.mpi / almanac.fpi;
  phys_diff = @(pmu)(phys_rat - fn_phys_rat(params, data, pmu));
  
  mu_max = max(data.mu);
  td.mu = fzero(phys_diff, [eps, mu_max]);
  
  td.meta.is_dummy = 1;
  td.meta.has_iso = data.meta.has_iso;
  td = calculate_predictions(params, td);
  
  res_fac = td.inf_fps / almanac.fpi;
  scale.a = res_fac * data.scale.a; % Inverse?
  scale.mu = td.mu / res_fac;
  scale.zp = fit_zp(params, data);

  mps_n_pred = @(pmu)(calc_mps_n(params, data, pmu));
  if data.meta.has_iso && (mps_n_pred(eps) < 0)
      scale.mu_min = fzero(mps_n_pred, [eps, mu_max]);
  else
      scale.mu_min = eps;
  end
  
  params.f0 = params.f0 / res_fac;
  params.B0 = params.B0 / res_fac;
  if data.meta.has_asq
      params.Dm = params.Dm / res_fac^2;
      params.Df = params.Df / res_fac^2;
  end
  if data.meta.has_iso
      params.zeta = params.zeta / res_fac^2;
  end
  if data.meta.needs_Dn
      params.Dn = params.Dn / res_fac^2;
  end
end

function fn_res = fn_phys_rat(params, data, mu)
  td.mu = mu;
  td.meta.is_dummy = 1;
  td.meta.has_iso = data.meta.has_iso;
  td = calculate_predictions(params, td);
  fn_res = sqrt(td.inf_mps2)/ td.inf_fps;
end

function mps_n = calc_mps_n(params, data, mu)
  td.mu = mu;
  td.meta.is_dummy = 1;
  td.meta.has_iso = data.meta.has_iso;
  td = calculate_predictions(params, td);
  mps_n = sqrt(td.inf_mps2_n);
end

function zp = fit_zp(params, data)
    if data.meta.num_betas == 1
        zp = data.meta.zp;
        return
    end
    
    zfac = ones(size(data.meta.zp));
    for idx = 1 : data.meta.num_betas - 1
        zfac(idx) = params.(data.meta.fn_zfac{idx});
    end
    zp = lscov(zfac, data.meta.zp, (data.meta.sd_zp).^-2);
end
