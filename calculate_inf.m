function data = calculate_inf(params, data)
    data.inf_mps2 = data.chimu_c .* (1 - data.xi_n .* data.lbar3_n);
    % Note that chimu_n contains lattice artifacts that should not be multiplying the logs!
    data.inf_mps2_n = data.chimu_c .* (1 - 2 .* data.xi_c .* data.lbar3_c + data.xi_n .* data.lbar3_n);
    if data.meta.needs_zeta
        data.inf_mps2_n = data.inf_mps2_n + params.zeta .* data.afac.^2;
    end
    if data.meta.has_iso
         data.inf_mps2_n = data.inf_mps2_n + params.zeta .* data.afac.^2 .* data.xi_n .* data.Xibar3_n;
    end

    data.inf_fps  = params.f0 * (1 + data.xi_c .* data.lbar4_c + data.xi_n .* data.lbar4_n);
    
    if data.meta.has_nnlo
        data.inf_mps2 = data.inf_mps2 + data.chimu_c .* ((17/2) .* data.xi_c.^2 .* ((1/51) .* (28 .* data.lbar1_c + 32 .* data.lbar2_c - 9 .* 28 .* data.lbar3_c + 49).^2)) + 4 .* data.xi_c.^2 .* params.km;
        data.inf_fps = data.inf_fps - params.f0  .* (5 .* data.xi_c.^2 .* ((1/30) .* (14 .* data.lbar1_c + 16 .* data.lbar2_c + 6 .* data.lbar3_c - 6 .* data.lbar4_c + 23).^2)) + 4 .* data.xi_c.^2 .* params.kf;
        if data.meta.has_iso
            data.inf_mps2_n = data.inf_mps2_n + data.chimu_c .* ((17/2) .* data.xi_c.^2 .* ((1/51) .* (28 .* data.lbar1_c + 32 .* data.lbar2_c - 9 .* 28 .* data.lbar3_c + 49).^2)) + 4 .* data.xi_c.^2 .* params.km;
        end
    end
end
