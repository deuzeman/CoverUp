function data = calculate_treelevel(params, data)
    data.chimu_c = 2 * params.B0 .* data.mu;
    data.chimu_n = data.chimu_c;
    if data.meta.has_iso
        data.chimu_n = data.chimu_n + params.zeta .* data.afac.^2;
    end
    data.xi_c = data.chimu_c ./ (16 * pi^2 * params.f0.^2);
    data.xi_n = data.chimu_n ./ (16 * pi^2 * params.f0.^2);
end
