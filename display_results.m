function display_results(params, data)
    global almanac;
    
    bcols = {'r', 'g', 'b'};
    
    [~, data] = chi(params, data);
    
    mu_max = ceil(5.0 * max(data.mu)) / 5.0;
    plot_data.mu = data.scale.mu : mu_max / 250 : mu_max;
    plot_data.meta.is_dummy = 1;
    plot_data.meta.has_iso  = data.meta.has_iso;
    plot_data.meta.needs_zeta = data.meta.needs_zeta;
  
    plot_data = calculate_predictions(params, plot_data);

    plot_data.inf_mps = sqrt(plot_data.inf_mps2);
    plot_data.inf_mps2_mu = plot_data.inf_mps2 ./ plot_data.mu;

    % Apply all corrections to the measured data
    data.fps = data.fps ./ data.fvol_fps - data.asq_fps;
    data.sd_fps = data.sd_fps ./ data.fvol_fps;
    
    data.mps2 = data.mps.^2 ./ data.fvol_mps2 - data.asq_mps2;
    data.sd_mps2 = 2.0 * data.sd_mps .* data.mps ./ data.fvol_mps2;
    
    data.mps2_n = data.mps_n.^2 ./ data.fvol_mps2_n - data.asq_mps2_n;
    data.sd_mps2_n = 2.0 * data.sd_mps_n .* data.mps_n ./ data.fvol_mps2_n;
    
    data.mps = sqrt(data.mps2);
    data.sd_mps = data.sd_mps ./ sqrt(data.fvol_mps2);
    
    data.mps2_mu = data.mps2 ./ data.mu;
    data.sd_mps2_mu = data.sd_mps2 ./ data.mu;
    
    figure('Name','Pion decay constant','NumberTitle','off');
    hold on;
    plot(plot_data.mu, plot_data.inf_fps, 'LineWidth', 2, 'Color', 'b');
    plot(data.scale.mu, almanac.fpi, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    for idx = 1 : data.meta.num_betas
        errorbar(data.mu(data.meta.indices{idx}), data.fps(data.meta.indices{idx}), ...
                 data.sd_fps(data.meta.indices{idx}), ...
                 'o', 'MarkerFaceColor', bcols{idx}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'Color', 'k');
    end
    hold off;
    
    figure('Name','Pion mass squared','NumberTitle','off');
    hold on;
    plot(plot_data.mu, plot_data.inf_mps2_mu, 'LineWidth', 2, 'Color', 'b');
    plot(data.scale.mu, almanac.mpi.^2 ./ data.scale.mu, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    for idx = 1 : data.meta.num_betas
        errorbar(data.mu(data.meta.indices{idx}), data.mps2_mu(data.meta.indices{idx}), ...
                 data.sd_mps2_mu(data.meta.indices{idx}), ...
                 'o', 'MarkerFaceColor', bcols{idx}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'Color', 'k');
    end
    hold off;
    
    figure('Name','Neutral pion mass squared','NumberTitle','off');
    hold on;
    plot(plot_data.mu, plot_data.inf_mps2_n, 'LineWidth', 2, 'Color', 'b');
    for idx = 1 : data.meta.num_betas
        errorbar(data.mu(data.meta.indices{idx}), data.mps2_n(data.meta.indices{idx}), ...
                 data.sd_mps2_mu(data.meta.indices{idx}), ...
                 'o', 'MarkerFaceColor', bcols{idx}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'Color', 'k');
    end
    hold off;
end
