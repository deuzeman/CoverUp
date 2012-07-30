function data = calculate_fvol(data)
    global opts;
    
    if data.meta.is_dummy
        return
    end
    
    lbar1 = repmat(data.lbar1_c, [1, 100]);
    lbar2 = repmat(data.lbar2_c, [1, 100]);
    lbar3 = repmat(data.lbar3_c, [1, 100]);
    lbar4 = repmat(data.lbar4_c, [1, 100]);

    mat_xi       = 2 * repmat(data.xi_c, 1, 100); % CDH convention for f0!
    lambda_c     = sqrt(data.chimu_c) .* data.L;
    sqn_lambda_c = lambda_c * sqrt(1:100); % NOTE Replaced by full pion mass, not LO! Probably wrong for GL correction scheme.
    prefactor    = repmat(multiplicity(1:100)', [length(data.mu), 1]) ./ sqn_lambda_c;

    switch upper(opts.fvol)
        case 'OFF'
            data.fvol_fps  = ones(size(data.mps));
            data.fvol_mps2 = ones(size(data.mps));
            data.fvol_mps2_n = data.fvol_mps2;
        case 'GL'
            data.fvol_fps  = 1 ./ (1 + 2 * data.xi_c .* delta1(lambda_c));
            data.fvol_mps2 = 1 ./ (1 - data.xi_c .* delta1(lambda_c));
            data.fvol_mps2_n = data.fvol_mps2;
        case 'BMW'
            mpl = data.mps .* data.L;
            data.fvol_fps  =  1 - 11.31937 * exp(-mpl) .* mpl.^(-1.5);
            data.fvol_mps2 = (1 + 7.706548 * exp(-mpl) .* mpl.^(-1.5)).^2;
            data.fvol_mps2_n = data.fvol_mps2;
        case 'CDH'
            I_mps_2 = -2 * besselk(1, sqn_lambda_c);
            I_mps_4 = (1 / 9) * (101 - 39 * pi + 72 * lbar1 + 48 * lbar2 - 45 * lbar3 - 36 * lbar4) .* besselk(1, sqn_lambda_c) ...
                      + (1 / 18) * (-476 + 183 * pi - 96 * lbar1 - 384 * lbar2) .* (besselk(2, sqn_lambda_c) ./ (sqn_lambda_c));

            I_fps_2 = -4 * besselk(1, sqn_lambda_c);
            I_fps_4 = ((1 / 36) * (58 - 87 * pi + 144 * lbar1 + 96 * lbar2 - 216 * lbar4) .*  besselk(1, sqn_lambda_c) ...
                      + (1 / 72) * (-2456 + 1173 * pi - 384 * lbar1 - 1536 * lbar2) .* (besselk(2, sqn_lambda_c) ./ (sqn_lambda_c)));

            data.fvol_fps  = 1 + sum(prefactor .* (mat_xi .* I_fps_2 + mat_xi.^2 .* I_fps_4), 2);
            data.fvol_mps2 = (1 - sum(0.5 * prefactor .* (mat_xi .* I_mps_2 + mat_xi.^2 .* I_mps_4), 2)).^2;
            data.fvol_mps2_n = data.fvol_mps2;
        case 'CWW'
            r = repmat(sqrt(data.chimu_n ./ data.chimu_c), 1, 100);

            sqn_lambda_n = sqn_lambda_c ./ r;

            Rcww_0 = Rcww(0, sqn_lambda_c);
            Rcww_1 = Rcww(1, sqn_lambda_c);
            Rcww_2 = Rcww(2, sqn_lambda_c);

            Rcww_0_n = Rcww(0, sqn_lambda_n, r);
            Rcww_1_n = Rcww(1, sqn_lambda_n, r);
            Rcww_2_n = Rcww(2, sqn_lambda_n, r);
            
            dRcww_0 = dRcww(0, sqn_lambda_c);
            dRcww_1 = dRcww(1, sqn_lambda_c);
            dRcww_2 = dRcww(2, sqn_lambda_c);

            dRcww_0_n = dRcww(0, sqn_lambda_n, r);
            dRcww_1_n = dRcww(1, sqn_lambda_n, r);
            dRcww_2_n = dRcww(2, sqn_lambda_n, r);

            B0 = 2 * besselk(1, sqn_lambda_c);
            B2 = 2 * besselk(2, sqn_lambda_c) ./ sqn_lambda_c;
            
            B0_n = 2 * r    .* besselk(1, sqn_lambda_n);
            B2_n = 2 * r.^3 .* besselk(2, sqn_lambda_n) ./ sqn_lambda_n;
            
            I_mps_2_c = 0;
            I_mps_2_n = -B0_n;

            I_mps_4_c = (( 8 / 3) * (lbar1 + lbar2) - 2 * lbar3 - (34 / 9)) .* B0 + ...
                        ((92 / 9) - (8 / 3) * lbar1 - 8 * lbar2) .* B2 + ...
                        + ( 1 / 3) * (11 * Rcww_0 + 20 * Rcww_1 - 32 * Rcww_2);
            I_mps_4_n = ((4 / 3) * lbar1 - 0.5 * lbar3 - 2 * lbar4 + (13 / 18)) .* B0_n + ...
                        ((20 / 9) - (8 / 3) * lbar2) .* B2_n + ...
                        + (2 / 3) * (Rcww_0_n + 2 * Rcww_1_n - 4 * Rcww_2_n);
                    
            I_fps_2_c = -B0;
            I_fps_2_n = -B0_n;
            
            I_fps_4_c = (-(8 / 9) + (4 / 3) * (lbar1 + lbar2) - 2 * lbar4) .* B0 + ...
                        ((92 / 9) - (8 / 3) * lbar1 - 8 * lbar2) .* B2 + ...
                        ( 1 / 3) * (2 * Rcww_0 - 8 * Rcww_1 - 32 * Rcww_2 + ...
                        - (11 / 2) * dRcww_0 + 10 * dRcww_1 + 16 * dRcww_2);
            I_fps_4_n = (( 1 / 9) + (2 / 3) * lbar1 - lbar4) .* B0_n + ...
                        ((20 / 9) - (8 / 3) * lbar2) .* B2_n + ...
                        ( 1 / 3) * (2 * Rcww_0_n + 4 * Rcww_1_n - 8 * Rcww_2_n ...
                         - dRcww_0_n - 2 * dRcww_1_n + 4 * dRcww_2_n);
            
            I_mps_pvci = -prefactor .* mat_xi.^3 .* Rpvci(sqn_lambda_n, r);
            I_fps_pvci = -0.75 * I_mps_pvci;

            data.fvol_fps    = 1 + sum(prefactor .* (mat_xi .* (I_fps_2_c + I_fps_2_n + mat_xi .* (I_fps_4_c + I_fps_4_n))) + I_fps_pvci, 2);
            data.fvol_mps2   = (1 - sum(0.5 * prefactor .* (mat_xi .* (I_mps_2_c + I_mps_2_n + mat_xi .* (I_mps_4_c + I_mps_4_n))) + I_mps_pvci, 2)).^2;
            
            I_mps_2_cn = 0;
            I_mps_n_2_c = 2 * I_mps_2_c;          
            I_mps_n_2_n = I_mps_2_cn - I_mps_2_n;
            
            I_mps_4_cn = (( 8 / 3) * (lbar1 + lbar2) - 2 * lbar3 - (34 / 9)) .* B0_n + ((92 / 9) - (8 / 3) * lbar1 - 8 * lbar2) .* B2_n + (1 / 3) * (11 * Rcww_0_n - 20 * Rcww_1_n - 32 * Rcww_2_n);
            I_mps_n_4_c = 2 * I_mps_4_c;
            I_mps_n_4_n = I_mps_4_cn - I_mps_4_n;
           
            data.fvol_mps2_n = (1 - sum(0.5 * prefactor .* (mat_xi .* (I_mps_n_2_c + I_mps_n_2_n + mat_xi .* (I_mps_n_4_c + I_mps_n_4_n))), 2)).^2;
        otherwise
            error(-1, ['Unsupported finite volume corrections: ' upper(opts.fvol)]);
    end
end

function res = multiplicity (n)
    mlock;
    persistent multn;
    if isempty (multn)
        multn = [   6, 12,  8,  6, 24,  24,  0, 12,  30,  24, 24,  8, 24,  48,  0,  6, 48,  36, 24, 24, ...
                   48, 24,  0, 24, 30,  72, 32,  0,  72,  48,  0, 12, 48,  48, 48, 30, 24,  72,  0, 24, ...
                   96, 48, 24, 24, 72,  48,  0,  8,  54,  84, 48, 24, 72,  96,  0, 48, 48,  24, 72,  0, ...
                   82, 96,  0,  6, 96,  96, 24, 48,  96,  48,  0, 36, 48, 120, 56, 24, 96,  48, 0,  24, ...
                  102, 48, 72, 48, 48, 120,  0, 24, 144, 120, 48,  0, 48,  96,  0, 24, 48, 108, 72, 30  ]';
    end
    res = multn (n);
end

function res = dRcww(k, sqn_lambda, r)
    if nargin == 2
        r = ones(size(sqn_lambda));
    end
    vec_sqn_lambda = reshape(sqn_lambda, [numel(sqn_lambda), 1]);
    vec_r = reshape(r, [numel(sqn_lambda), 1]);
    fun = @(y)dintv(k, vec_sqn_lambda, vec_r, y);
    if mod(k, 2)
        res = reshape(imag(quadv(fun, -50, 50)), size(sqn_lambda));
    else
        res = reshape(real(quadv(fun, -50, 50)), size(sqn_lambda));
    end  
end

function res = Rcww(k, sqn_lambda, r)
    if nargin == 2
        r = ones(size(sqn_lambda));
    end
    vec_sqn_lambda = reshape(sqn_lambda, [numel(sqn_lambda), 1]);
    vec_r = reshape(r, [numel(sqn_lambda), 1]);
    fun = @(y)intv(k, vec_sqn_lambda, vec_r, y);
    if mod(k, 2)
        res = reshape(imag(quadv(fun, -50, 50)), size(sqn_lambda));
    else
        res = reshape(real(quadv(fun, -50, 50)), size(sqn_lambda));
    end   
end

function res = Rpvci(sqn_lambda, r)
    if nargin == 2
        r = ones(size(sqn_lambda));
    end

    w = 0.5 * r.^-1;
    vec_sqn_lambda = reshape(sqn_lambda, [numel(sqn_lambda), 1]);
    vec_w = reshape(w, [numel(sqn_lambda), 1]);
    
    fun = @(y)(exp(-vec_sqn_lambda .* sqrt(1 + y.^2)) ./ (y.^2 + vec_w.^2));
    res = reshape(quadv(fun, -50, 50), size(sqn_lambda));
    
    res = pi * exp(-sqn_lambda .* sqrt(1 - w.^2)) - 0.25 * res;
end

function res = intv(k, vec_sqn_lambda, vec_r, y)
    res = y.^k .* exp(-sqrt(1 + (y ./ vec_r).^2) .* vec_sqn_lambda) .* g(2 + 2i * y);
end

function res = dintv(k, vec_sqn_lambda, vec_r, y)
    res = y.^k .* exp(-sqrt(1 + (y ./ vec_r).^2) .* vec_sqn_lambda) .* dg(2 + 2i * y);
end

function res = g(x)
    res = 2 + sqrt((x - 4) ./ x) .* log(1 + 0.5 * sqrt(x - 4) .* sqrt(x) - 0.5 * x);
end

function res = dg(x)
    res = 4 - x + 2 * sqrt((x - 4) ./ x) .* log(1 + 0.5 * sqrt(x - 4) .* sqrt(x) - 0.5 * x);
    res = res ./ ((x - 4) .* x);
end

function res = delta1(sqn_lambda)
    % shape function as defined in hep-lat/0804.0473v1 equation C5 on page 84
    multn = repmat(multiplicity(1 : 100)', size(sqn_lambda, 1), 1);
    res = sum((4 * multn ./ sqn_lambda) .* besselk(1, sqn_lambda), 2);
end

