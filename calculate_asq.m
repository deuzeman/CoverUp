function data = calculate_asq(params, data)   
    if data.meta.is_dummy
        return
    end
    
    data.asq_mps2   = zeros(size(data.mu));
    data.asq_fps    = zeros(size(data.mu));
    data.asq_mps2_n = zeros(size(data.mu));

    if data.meta.has_asq 
        data.asq_mps2   = data.chimu_c .* params.Dm .* data.afac.^2;
        data.asq_fps    = params.f0 .* params.Df .* data.afac.^2;
    end
    if data.meta.needs_Dn
        data.asq_mps2_n = params.Dn .* data.afac.^4;
    end

end
