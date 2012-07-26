function data = calculate_predictions(params, data)
    data = rescale_data(params, data);
    data = calculate_treelevel(params, data);
    data = renormalize_lecs(params, data);
    data = calculate_inf(params, data);
    data = calculate_fvol(params, data);
    data = calculate_asq(params, data);    
end
