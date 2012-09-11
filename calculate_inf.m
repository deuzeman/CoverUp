function data = calculate_inf(data)
    data.inf_mps2 = data.chimu_c .* (1 - data.xi_n .* data.lbar3_n);
    data.inf_fps  = data.params.f0 * (1 + data.xi_c .* data.lbar4_c + data.xi_n .* data.lbar4_n);
    
    if data.meta.is_dummy
        data.inf_mps2_n = data.inf_mps2;
        data.inf_mps = sqrt(data.inf_mps2);
        data.inf_mps_n = sqrt(data.inf_mps2_n);
        return
    end

    % Note that chimu_n contains lattice artifacts that should not be multiplying the logs!
    data.inf_mps2_n = data.chimu_c .* (1 - 2 .* data.xi_c .* data.lbar3_c + data.xi_n .* data.lbar3_n);
    
    if data.meta.needs_zeta
        data.inf_mps2_n = data.inf_mps2_n - data.params.zeta .* data.afac.^2;
    end
    if data.meta.has_iso
        data.inf_mps2_n = data.inf_mps2_n - data.params.zeta .* data.xi_n .* data.Xibar3_n .* data.afac.^2;
    end
    
    if data.meta.has_nnlo
        data.inf_mps2 = data.inf_mps2 + data.chimu_c .* ((17/2) .* data.xi_c.^2 .* ((1/51) .* (28 .* data.lbar1_c + 32 .* data.lbar2_c - 9 .* 28 .* data.lbar3_c + 49).^2)) + 4 .* data.xi_c.^2 .* data.params.km;
        data.inf_fps = data.inf_fps - data.params.f0  .* (5 .* data.xi_c.^2 .* ((1/30) .* (14 .* data.lbar1_c + 16 .* data.lbar2_c + 6 .* data.lbar3_c - 6 .* data.lbar4_c + 23).^2)) + 4 .* data.xi_c.^2 .* data.params.kf;
        if data.meta.has_iso
            data.inf_mps2_n = data.inf_mps2_n + data.chimu_c .* ((17/2) .* data.xi_c.^2 .* ((1/51) .* (28 .* data.lbar1_c + 32 .* data.lbar2_c - 9 .* 28 .* data.lbar3_c + 49).^2)) + 4 .* data.xi_c.^2 .* data.params.km;
        end
    end
    
    % Protect against negative neutral pion mass by flattening
    data.inf_mps2_n = max(1e-3, data.inf_mps2_n);
    
    data.inf_mps = sqrt(data.inf_mps2);
    data.inf_mps_n = sqrt(data.inf_mps2_n);
end
