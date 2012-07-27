function data = renormalize_lecs(params, data)
    global almanac;
    
    % If renormalization should happen in the continuum limit, we have only one value!
    
    data.lbar1_c = params.l1 + 2 * log(almanac.mpi) - log(data.chimu_c);
    data.lbar1_n = params.l1 + 2 * log(almanac.mpi) - log(data.chimu_c);

    data.lbar2_c = params.l2 + 2 * log(almanac.mpi) - log(data.chimu_c);
    data.lbar2_n = params.l2 + 2 * log(almanac.mpi) - log(data.chimu_c);

    data.lbar3_c = params.l3 + 2 * log(almanac.mpi) - log(data.chimu_c);
    data.lbar3_n = params.l3 + 2 * log(almanac.mpi) - log(data.chimu_c);

    data.lbar4_c = params.l4 + 2 * log(almanac.mpi) - log(data.chimu_c);
    data.lbar4_n = params.l4 + 2 * log(almanac.mpi) - log(data.chimu_c);
    
    if data.meta.has_iso
        data.Xibar3_n = params.Xi3 + 2 * log(almanac.mpi) - log(data.chimu_c);
    end
end
