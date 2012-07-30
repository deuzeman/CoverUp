function display_results(params, data)
    global almanac;
    global opts;
    
    bcols = {'r', 'g', 'b'};
    
    [~, data] = chi(params, data);
    
    mu_max = ceil(5.0 * max(data.mu)) / 5.0;
    plot_data.mu = data.scale.mu : mu_max / 250 : mu_max;
    plot_data.meta.is_dummy = 1;
    plot_data.meta.has_iso  = data.meta.has_iso;
    plot_data.meta.needs_zeta = data.meta.needs_zeta;
    plot_data.meta.has_nnlo = data.meta.has_nnlo;
  
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
    
    figure('Name','Pion mass','NumberTitle','off');
    hold on;
    plot(plot_data.mu, sqrt(plot_data.inf_mps2), 'LineWidth', 2, 'Color', 'b');
    plot(data.scale.mu, almanac.mpi, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    for idx = 1 : data.meta.num_betas
        errorbar(data.mu(data.meta.indices{idx}), data.mps(data.meta.indices{idx}), ...
                 data.sd_mps(data.meta.indices{idx}), ...
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
    
    if data.meta.has_iso
        figure('Name','Neutral pion mass','NumberTitle','off');
        hold on;
        plot(plot_data.mu, sqrt(plot_data.inf_mps2_n), 'LineWidth', 2, 'Color', 'b');
        for idx = 1 : data.meta.num_betas
            errorbar(data.mu(data.meta.indices{idx}), sqrt((data.mps2_n(data.meta.indices{idx}) + params.zeta * data.afac(data.meta.indices{idx}).^2 - data.asq_mps2_n(data.meta.indices{idx})) ./ data.fvol_mps2_n(data.meta.indices{idx})), ...
                     data.sd_mps_n(data.meta.indices{idx}), ...
                     'o', 'MarkerFaceColor', bcols{idx}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'Color', 'k');
        end
        hold off;
    end
    fprintf('Fitting report\n');
    fprintf('====================================\n');
    fprintf('  Settings\n');
    fprintf('    O(a^2) artefacts turned %s.\n', opts.asq);
    fprintf('    Finite volume effects calculated through %s.\n', opts.fvol);
    fprintf('    Isospin splitting effects turned %s.\n', opts.iso);
    fprintf('    Prior inclusion set to %s.\n', opts.priors);
    fprintf('    NNLO terms set to %s.\n', opts.nnlo);
    fprintf('====================================\n');
    fprintf('  Results\n');
    fprintf('    a                 : %7.4f fm\n', data.scale.a * almanac.mev_fm);
    fprintf('    f_0               : %7.3f MeV\n', params.f0);
    fprintf('    f_pi / f_0        : %7.3f\n\n', almanac.fpi / params.f0);
    fprintf('    B_0               : %7.3f MeV\n', params.B0);
    fprintf('    2 B_0 mu / m_pi^2 : %7.3f\n\n', 2 * params.B0 * data.scale.mu / almanac.mpi^2);
    fprintf('    \\bar{l}_1         : %7.3f\n', params.l1);
    fprintf('    \\bar{l}_2         : %7.3f\n', params.l2);
    fprintf('    \\bar{l}_3         : %7.3f\n', params.l3);
    fprintf('    \\bar{l}_4         : %7.3f\n\n', params.l4);
    if data.meta.has_asq
        fprintf('    D_m               : %7.3f fm^2\n', params.Dm / (data.scale.a * almanac.mev_fm)^2);
        fprintf('    D_f               : %7.3f fm^2\n', params.Df / (data.scale.a * almanac.mev_fm)^2);
        if data.meta.needs_Dn
            fprintf('    D_n              : %7.3f fm^2\n', params.Dn / (data.scale.a * almanac.mev_fm)^2);
        end
    end

    if data.meta.needs_zeta
        fprintf('    \\zeta             : %7.3e fm^2\n', params.zeta / (data.scale.a * almanac.mev_fm)^2);
    end
    if data.meta.has_iso
        fprintf('    \\Xi_3             : %7.3f\n', params.Xi3);
    end
    fprintf('====================================\n');
    fprintf('  Details\n');
    fprintf('    Total chi squared  : %5.2f\n', sum(s2v(data.chi).^2));
    fprintf('      From pion mass   : %5.2f\n', sum((data.chi.mps).^2));
    fprintf('      From pion decay  : %5.2f\n', sum((data.chi.fps).^2));
    if data.meta.has_iso
        fprintf('      From neutral pion : %5.2f\n', sum((data.chi.mps2_n).^2));
    end
    if ~strcmpi(opts.priors, 'OFF')
        fprintf('      From priors      : %5.2f\n', sum((data.chi.priors).^2));
    end
    fprintf('    Chi squared / dof  : %5.2f\n', sum(s2v(data.chi).^2) / (length(s2v(data.chi)) - 6 - data.meta.has_asq * 2 - data.meta.needs_Dn - data.meta.needs_zeta - data.meta.has_iso));
end

function vec = s2v(st)
    f = fieldnames(st);
    vec = [];
    for ctr = 1 : length(f)
        vec = [vec; st.(f{ctr})]; %#ok
    end
end
