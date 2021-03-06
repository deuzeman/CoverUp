function data = calculate_treelevel(data)   
    data.chimu_c = 2 * data.params.B0 .* data.mu;
    data.chimu_n = data.chimu_c;
    data.xi_c = data.chimu_c ./ (16 * pi^2 * data.params.f0.^2);
    
    if data.meta.is_dummy
        data.xi_n = data.xi_c;
        return
    end
   
    if data.meta.needs_zeta
        data.chimu_n = data.chimu_n + 2.0 * data.params.zeta .* data.afac.^2;
    end

    data.xi_n = data.chimu_n ./ (16 * pi^2 * data.params.f0.^2);
    data.xi_n = max(1e-3, data.xi_n);
end
