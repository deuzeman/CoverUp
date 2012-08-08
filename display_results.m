function data = display_results(data)
    global almanac;
    global opts;
    
    bcols = {'r', 'g', 'b'};

    % Calculate clean reference data (a = 0, L = infinite)
    clean_data = data;
    clean_data.meta.clean = 1;
    clean_data = calculate_predictions(clean_data);   
    
    [~, data] = chi(data.params, data);
    
    % Calculate lines for plotting
    mu_max = ceil(5.0 * max(data.mu)) / 5.0;
    plot_data.mu = eps : mu_max / 250 : mu_max;
    plot_data.meta = data.meta;
    plot_data.meta.clean = 1;
    plot_data.params = data.params;
    
    plot_data = calculate_predictions(plot_data);

    plot_data.inf_mps = sqrt(plot_data.inf_mps2);
    plot_data.inf_mps2_mu = plot_data.inf_mps2 ./ plot_data.mu;
        
    % Calculate deviations
    data.fps_corr = clean_data.inf_fps + data.dev.fps;
    data.mps_corr = clean_data.inf_mps + data.dev.mps;
    data.mps_n_corr = clean_data.inf_mps_n + data.dev.mps_n;
    
    figure('Name', sprintf('Fit %d: Pion decay constant', opts.fit_cnt), 'NumberTitle','off');
    hold on;
    plot(plot_data.mu, plot_data.inf_fps, 'LineWidth', 2, 'Color', 'b');
    plot(data.scale.mu, almanac.fpi, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    for idx = 1 : data.meta.num_betas
        errorbar(data.mu(data.meta.indices{idx}), data.fps_corr(data.meta.indices{idx}), ...
                 data.sd_fps(data.meta.indices{idx}), ...
                 'o', 'MarkerFaceColor', bcols{idx}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'Color', 'k');
    end
    hold off;
    
    figure('Name', sprintf('Fit %d: Pion mass', opts.fit_cnt), 'NumberTitle','off');
    hold on;
    plot(plot_data.mu, sqrt(plot_data.inf_mps2), 'LineWidth', 2, 'Color', 'b');
    plot(data.scale.mu, almanac.mpi, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    for idx = 1 : data.meta.num_betas
        errorbar(data.mu(data.meta.indices{idx}), data.mps_corr(data.meta.indices{idx}), ...
                 data.sd_mps(data.meta.indices{idx}), ...
                 'o', 'MarkerFaceColor', bcols{idx}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'Color', 'k');
    end
    hold off;
    
    figure('Name', sprintf('Fit %d: Pion mass squared over mu', opts.fit_cnt), 'NumberTitle','off');
    hold on;
    plot(plot_data.mu, plot_data.inf_mps2 ./ plot_data.mu, 'LineWidth', 2, 'Color', 'b');
    plot(data.scale.mu, almanac.mpi.^2 ./ data.scale.mu, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
    for idx = 1 : data.meta.num_betas
        errorbar(data.mu(data.meta.indices{idx}), data.mps_corr(data.meta.indices{idx}).^2 ./ data.mu(data.meta.indices{idx}), ...
                 2 * data.mps(data.meta.indices{idx}) .* data.sd_mps(data.meta.indices{idx}) ./ data.mu(data.meta.indices{idx}), ...
                 'o', 'MarkerFaceColor', bcols{idx}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'Color', 'k');
    end
    hold off;
    
    if data.meta.has_iso
        figure('Name', sprintf('Fit %d: Neutral pion mass', opts.fit_cnt), 'NumberTitle','off');
        hold on;
        plot(plot_data.mu, sqrt(plot_data.inf_mps2_n), 'LineWidth', 2, 'Color', 'b');
        for idx = 1 : data.meta.num_betas
            errorbar(data.mu(data.meta.indices{idx}), data.mps_n_corr(data.meta.indices{idx}), ...
                     data.sd_mps_n(data.meta.indices{idx}), ...
                     'o', 'MarkerFaceColor', bcols{idx}, 'MarkerEdgeColor', 'k', 'LineWidth', 2, 'Color', 'k');
        end
        hold off;
    end
    fprintf('Fitting report for fit %d\n', opts.fit_cnt);
    fprintf('=========================================================\n');
    fprintf('  Settings\n');
    fprintf('    O(a^2) artefacts turned %s.\n', opts.asq);
    fprintf('    Finite volume effects calculated through %s.\n', opts.fvol);
    fprintf('    Isospin splitting effects turned %s.\n', opts.iso);
    fprintf('    Prior inclusion set to %s.\n', opts.priors);
    fprintf('    NNLO terms set to %s.\n', opts.nnlo);
    fprintf('    Run with %d bootstrap samples.\n', opts.nboot);
    fprintf('    Data writeout is turned %s.\n', opts.write);
    fprintf('=========================================================\n');
    fprintf('  Results\n');
    for ctr = 1 : data.meta.num_betas - 1
        fprintf('    a_{%4.2f}          : %7.4f +/- %7.4f fm\n', data.meta.betas(ctr), data.scale.a * almanac.mev_fm * data.params.(data.meta.fn_afac{ctr}), ...
                                                                                       data.sd_scale.a * almanac.mev_fm * data.params.(data.meta.fn_afac{ctr}));
    end
    fprintf('    a_{%4.2f}          : %7.4f +/- %7.4f fm\n', data.meta.betas(end), data.scale.a * almanac.mev_fm, data.sd_scale.a * almanac.mev_fm);
    fprintf('    f_0               : %7.3f +/- %7.3f MeV\n', data.params.f0, data.sd_params.f0);
    fprintf('    f_pi / f_0        : %7.3f +/- %7.3f\n\n', almanac.fpi / data.params.f0, almanac.fpi * data.sd_params.f0 / data.params.f0^2);
    fprintf('    B_0               : %7.2f +/- %7.2f MeV\n', data.params.B0, data.sd_params.B0);
    fprintf('    \\mu / Z_p         : %7.3f +/- %7.3f MeV\n', data.scale.mu, data.sd_scale.mu);
    fprintf('    2 B_0 mu / m_pi^2 : %7.3f +/- %7.3f\n\n', 2 * data.params.B0 * data.scale.mu / almanac.mpi^2, 2 * data.sd_params.B0 * data.scale.mu / almanac.mpi^2);
    fprintf('    \\bar{l}_1         : %7.3f +/- %7.3f\n', data.params.l1, data.sd_params.l1);
    fprintf('    \\bar{l}_2         : %7.3f +/- %7.3f\n', data.params.l2, data.sd_params.l2);
    fprintf('    \\bar{l}_3         : %7.3f +/- %7.3f\n', data.params.l3, data.sd_params.l3);
    fprintf('    \\bar{l}_4         : %7.3f +/- %7.3f\n\n', data.params.l4, data.sd_params.l4);
    if data.meta.has_asq
        fprintf('    D_m                : %7.3f +/- %7.3f fm^{-2}\n', data.params.Dm / (data.scale.a * almanac.mev_fm)^2, data.sd_params.Dm / (data.scale.a * almanac.mev_fm)^2);
        fprintf('    D_f                : %7.3f +/- %7.3f fm^{-2}\n', data.params.Df / (data.scale.a * almanac.mev_fm)^2, data.sd_params.Df / (data.scale.a * almanac.mev_fm)^2);
        if data.meta.needs_Dn
            fprintf('    D_n               : %7.3f +/- %7.3f fm^{-2}\n', data.params.Dn / (data.scale.a * almanac.mev_fm)^2, data.sd_params.Dn / (data.scale.a * almanac.mev_fm)^2);
        end
    end

    if data.meta.needs_zeta
        fprintf('    \\zeta             : %7.3f +/- %7.3f GeV^2 fm^{-2}\n', data.params.zeta / (1000 * data.scale.a * almanac.mev_fm)^2, data.sd_params.zeta / (1000 * data.scale.a * almanac.mev_fm)^2);
    end
    if data.meta.has_iso
        fprintf('    \\Xi_3             : %7.3f +/- %7.3f\n', data.params.Xi3, data.sd_params.Xi3);
    end
    fprintf('=========================================================\n');
    fprintf('  Fit quality details\n');
    fprintf('    From pion mass    : %7.2f\n', sum((data.chi.mps).^2));
    fprintf('    From pion decay   : %7.2f\n', sum((data.chi.fps).^2));
    if data.meta.has_iso
        fprintf('    From neutral pion : %7.2f\n', sum((data.chi.mps_n).^2));
    end
    if ~strcmpi(opts.priors, 'OFF')
        fprintf('    From non-pion     : %7.2f\n', sum((data.chi.priors).^2));
    end
    
    chisq   = sum(s2v(data.chi).^2);
    dof_chi = length(s2v(data.chi)); 
    dof_par = 6 + data.meta.has_asq * 2 + data.meta.needs_Dn + data.meta.needs_zeta + data.meta.has_iso;
    fprintf('    ---------------------------\n');
    fprintf('    Chi squared / dof : %7.2f / (%d - %d) = %7.2f\n', chisq, dof_chi, dof_par, chisq / (dof_chi - dof_par));
    fprintf('=========================================================\n');
    % Prepare data for writeout
    if strcmpi(opts.write, 'ON')
        data_write = [data.raw.beta, data.mu, data.afac .* data.scale.a, data.L, ...
                      data.mps, data.sd_mps, data.fps, data.sd_fps, data.mps_n, data.sd_mps_n, ...
                      data.fvol_mps2, data.asq_mps2, data.fvol_fps, data.asq_fps, data.fvol_mps2_n, data.asq_mps2_n, ...
                      data.mn,  data.sd_mn,  data.mk, data.sd_mk, data.md, data.sd_md];
        plot_write = [plot_data.mu', plot_data.inf_mps', plot_data.inf_fps', sqrt(plot_data.inf_mps2_n)'];

        filename_data = sprintf('Results/fit_%d.dat', opts.fit_cnt);
        filename_plot = sprintf('Results/fit_%d_plot.dat', opts.fit_cnt);

        file_data = fopen(filename_data, 'w');
        fprintf(file_data, '\tbeta\tmu\ta\tL\tmps\tsd_mps\tfps\tsd_fps\tmps_n\tsd_mps_n\tfvol_mps2\tasq_mps2\tfvol_fps\tasq_fps\tfvol_mps2_n\tasq_mps2_n\tmn\tsd_mn\tmk\tsd_mk\tmd\tsd_md\n');
        for idx = 1 : size(data_write, 1)
            fprintf(file_data, '%s\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\t%E\n', data.raw.ensemble(idx,:), data_write(idx, :));
        end
        fclose(file_data);

        file_plot = fopen(filename_plot,'w');
        fprintf(file_plot, 'mu\tmps\tfps\tmps_n\n');
        for idx = 1 : size(plot_write, 1)
            fprintf(file_plot, '%E\t%E\t%E\t%E\n', plot_write(idx, :));
        end
        fclose(file_plot);
    end
    
    opts.fit_cnt = opts.fit_cnt + 1;
end

function vec = s2v(st)
    f = fieldnames(st);
    vec = [];
    for ctr = 1 : length(f)
        vec = [vec; st.(f{ctr})]; %#ok
    end
end
