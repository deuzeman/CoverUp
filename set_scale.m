function data = set_scale(data)  
  global almanac;
  
  phys_rat = almanac.mpi / almanac.fpi;
  phys_diff = @(pmu)(real(phys_rat - fn_phys_rat(data, pmu)));
  
  mu_max = max(data.mu);
  td.mu = fzero(phys_diff, [eps, mu_max]);
  
  td.meta = data.meta;
  td.meta.is_dummy = 1;
  td.params = data.params;
  td = calculate_predictions(td);
  
  res_fac = td.inf_fps / almanac.fpi;
  data.scale.a = res_fac * data.scale.a;
  data.scale.mu = td.mu / res_fac;
  data = fit_zp(data);

  mps_n_pred = @(pmu)(calc_mps_n(data, pmu));
  if data.meta.has_iso && (mps_n_pred(eps) < 0)
      data.scale.mu_min = fzero(mps_n_pred, [eps, mu_max]);
  else
      data.scale.mu_min = eps;
  end
  
  data.params.f0 = data.params.f0 / res_fac;
  data.params.B0 = data.params.B0 / res_fac;
  if data.meta.has_asq
      data.params.Dm = data.params.Dm / res_fac^2;
      data.params.Df = data.params.Df / res_fac^2;
  end
  if data.meta.has_iso
      data.params.zeta = data.params.zeta / res_fac^2;
  end
  if data.meta.needs_Dn
      data.params.Dn = data.params.Dn / res_fac^4;
  end
end

function fn_res = fn_phys_rat(data, mu)
  td.mu = mu;
  td.meta = data.meta;
  td.meta.is_dummy = 1;
  td.params = data.params;
  td = calculate_predictions(td);
  fn_res = td.inf_mps / td.inf_fps;
end

function mps_n = calc_mps_n(data, mu)
  td.mu = mu;
  td.meta = data.meta;
  td.meta.is_dummy = 1;
  td.params = data.params;
  td = calculate_predictions(td);
  mps_n = td.inf_mps_n;
end

function data = fit_zp(data)
    if data.meta.num_betas == 1
        data.scale.zp = data.meta.zp;
        return
    end
    
    zfac = ones(size(data.meta.zp));
    for idx = 1 : data.meta.num_betas - 1
        zfac(idx) = data.params.(data.meta.fn_zfac{idx});
    end
    data.scale.zp = lscov(zfac, data.meta.zp, (data.meta.sd_zp).^-2);
end
