function data = renormalize_lecs(data)
    global almanac;
      
    if data.meta.has_l12
        data.lbar1_c = data.params.l1 + 2 * log(almanac.mpi) - log(data.chimu_c);
        data.lbar1_n = data.params.l1 + 2 * log(almanac.mpi) - log(data.chimu_n);

        data.lbar2_c = data.params.l2 + 2 * log(almanac.mpi) - log(data.chimu_c);
        data.lbar2_n = data.params.l2 + 2 * log(almanac.mpi) - log(data.chimu_n);
    end

    data.lbar3_c = data.params.l3 + 2 * log(almanac.mpi) - log(data.chimu_c);
    data.lbar3_n = data.params.l3 + 2 * log(almanac.mpi) - log(data.chimu_n);

    data.lbar4_c = data.params.l4 + 2 * log(almanac.mpi) - log(data.chimu_c);
    data.lbar4_n = data.params.l4 + 2 * log(almanac.mpi) - log(data.chimu_n);
    
    if data.meta.has_iso
        data.Xibar3_n = data.params.Xi3 + 2 * log(almanac.mpi) - log(data.chimu_n);
    end
end
