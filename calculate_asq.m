function data = calculate_asq(data)

    data.asq_mps2   = zeros(size(data.mu));
    data.asq_fps    = zeros(size(data.mu));
    data.asq_mps2_n = zeros(size(data.mu));
    
    if data.meta.is_dummy
        return;
    end

    if data.meta.has_asq 
        data.asq_mps2   = data.chimu_c .* data.params.Dm .* data.afac.^2;
        data.asq_fps    = data.params.f0 .* data.params.Df .* data.afac.^2;
    end
    if data.meta.needs_Dn
        data.asq_mps2_n = data.params.zeta .* data.params.Dn .* data.afac.^4;
    end

end
