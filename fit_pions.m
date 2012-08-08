function data =  fit_pions(data)
    % Using the estimated scales, set up a manual iteration
    for sc_iter = 1 : 20
        % Merge in the scale
        data = multiply_a(data);
        
        % Perform the fit
        data = lsq(data);
        
        % Recalibrate the scale
        old_scale = data.scale;
        data = set_scale(data);
    
        if sqrt(((data.scale.a - old_scale.a)^2 + (data.scale.zp - old_scale.zp)^2)) < 1e-8
            % Note that one would need to retune the scales of the data if this convergence criterion is set coarsely!
            break;
        end
    end
end

function data = multiply_a(data)
    data.mu  = data.raw.a_mu / (data.scale.a * data.scale.zp);               
    data.L = data.raw.L_a * data.scale.a;

    data.fps    = data.raw.a_fps / data.scale.a;
    data.sd_fps = data.raw.sd_a_fps / data.scale.a;
    data.mps    = data.raw.a_mps / data.scale.a;
    data.sd_mps = data.raw.sd_a_mps / data.scale.a;
    data.mps_n    = data.raw.a_mps_n / data.scale.a;
    data.sd_mps_n = data.raw.sd_a_mps_n / data.scale.a;
    
    data.mn       = data.raw.a_mn / data.scale.a;
    data.sd_mn    = data.raw.sd_a_mn / data.scale.a;
    data.mk       = data.raw.a_mk / data.scale.a;
    data.sd_mk    = data.raw.sd_a_mk / data.scale.a;
    data.md       = data.raw.a_md / data.scale.a;
    data.sd_md    = data.raw.sd_a_md / data.scale.a;
end
