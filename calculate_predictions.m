function data = calculate_predictions(data)
    data = rescale_data(data);
    data = calculate_treelevel(data);
    data = renormalize_lecs(data);
    data = calculate_inf(data);
    data = calculate_fvol(data);
    data = calculate_asq(data);    
end
