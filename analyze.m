function data = analyze(data)
    global opts;

    if strcmpi(opts.wipe, 'ON')
        close all;
    end
    
    data = prepare_data(data);

    data = fitpions(data);
    
    display_results(data);
end
