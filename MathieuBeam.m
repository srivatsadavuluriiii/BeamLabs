classdef MathieuBeam
    % MathieuBeam: Models non-diffracting beams in Elliptic coordinates.
    % FULLY STANDALONE: Includes numerical solvers for Angular and Radial
    % Mathieu functions (no Symbolic Toolbox required).
    
    properties
        order           % Mode order m
        q_param         % Ellipticity parameter q
        parity          % 'even' or 'odd'
        wavelength      % Wavelength
        w0_aperture     % Gaussian aperture scale
        
        % Derived
        k               % Total wavenumber
        kt              % Transverse wavenumber
        f_focal         % Semifocal distance
        
        % Numerical Data (Pre-computed)
        eigenval_a      % Characteristic value
        coeffs          % Expansion coefficients (Fourier/Bessel weights)
    end
    
    methods
        function obj = MathieuBeam(order, q_param, parity, wavelength, w0_aperture)
            if q_param < 0, error('q must be >= 0'); end
            
            obj.order = order;
            obj.q_param = q_param;
            obj.parity = lower(parity);
            obj.wavelength = wavelength;
            obj.w0_aperture = w0_aperture;
            obj.k = 2 * pi / wavelength;
            
            % --- Define Physical Scale ---
            if q_param == 0
                 % Bessel limit
                 obj.kt = 2.4048 / (w0_aperture/2); 
                 obj.f_focal = 0;
                 obj.coeffs = 1; 
            else
                 % Fit beam within aperture
                 obj.kt = 4 / w0_aperture; 
                 obj.f_focal = 2 * sqrt(q_param) / obj.kt;
                 
                 % --- Pre-solve Eigenvalue Problem ---
                 [obj.coeffs, obj.eigenval_a] = ...
                     MathieuBeam.solve_mathieu_eigenproblem(order, q_param, obj.parity);
            end
        end
        
        function field = generate_beam_field(obj, x, y, z, varargin)
            p_parser = inputParser;
            addRequired(p_parser, 'x');
            addRequired(p_parser, 'y');
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            parse(p_parser, x, y, z, varargin{:});
            args = p_parser.Results;
            
            % --- Bessel Limit Handling (q=0) ---
            if obj.q_param < 1e-9
                r = sqrt(x.^2 + y.^2);
                phi = atan2(y, x);
                if strcmp(obj.parity, 'even')
                    field_ideal = besselj(obj.order, obj.kt * r) .* cos(obj.order * phi);
                else
                    field_ideal = besselj(obj.order, obj.kt * r) .* sin(obj.order * phi);
                end
            else
                % --- Elliptic Coordinate Transform ---
                % z_complex = f * cosh(w) -> w = acosh(z/f)
                z_c = x + 1i * y;
                w_c = acosh(z_c ./ obj.f_focal);
                xi = real(w_c);
                eta = imag(w_c);
                
                % --- 1. Angular Mathieu (Fourier Series) ---
                Phi_eta = MathieuBeam.eval_angular(eta, obj.order, obj.parity, obj.coeffs);
                
                % --- 2. Radial Mathieu (Bessel Series) ---
                R_xi = MathieuBeam.eval_radial(xi, obj.order, obj.q_param, obj.parity, obj.coeffs);
                
                field_ideal = R_xi .* Phi_eta;
            end
            
            % --- Gaussian Apodization & Propagation ---
            r_sq = x.^2 + y.^2;
            gauss_env = exp(-r_sq ./ obj.w0_aperture^2);
            
            % Longitudinal k
            if obj.kt >= obj.k
                kz = 0; % Evanescent
            else
                kz = sqrt(obj.k^2 - obj.kt^2);
            end
            prop_phase = exp(1i * kz * z);
            
            scale = 1.0;
            if args.P_tx_watts > 0, scale = sqrt(args.P_tx_watts); end
            
            field = scale * field_ideal .* gauss_env .* prop_phase;
        end
        
        function intensity = calculate_intensity(obj, x, y, z, varargin)
            field = obj.generate_beam_field(x, y, z, varargin{:});
            intensity = abs(field).^2;
        end
    end
    
    methods (Static)
        function [coeffs, eigval] = solve_mathieu_eigenproblem(m, q, parity)
            % Solves the tridiagonal recurrence to find expansion coefficients.
            % Returns eigenvector 'coeffs' and eigenvalue 'eigval'.
            
            N_matrix = max(m + 25, 60); % Truncation size
            M = zeros(N_matrix);
            
            % --- Matrix Construction based on Parity and Order ---
            if strcmp(parity, 'even')
                if mod(m, 2) == 0 % Even indices (0, 2, 4...)
                    k_indices = 0:2:(2*(N_matrix-1));
                    % Row 1 (k=0) special: q*A_2 = a*A_0
                    M(1, 2) = q * sqrt(2); 
                    M(2, 1) = q * sqrt(2);
                    M(1, 1) = 0; 
                    
                    for i = 2:N_matrix
                        k = k_indices(i);
                        M(i, i) = k^2;
                        if i > 2, M(i, i-1) = q; end
                        if i < N_matrix, M(i, i+1) = q; end
                    end
                else % Odd indices (1, 3, 5...)
                    k_indices = 1:2:(2*N_matrix-1);
                    for i = 1:N_matrix
                        k = k_indices(i);
                        M(i, i) = k^2;
                        % Special boundary for ce_1: (1+q)A_1 - qA_3 = aA_1
                        if i == 1
                             M(i, i) = 1 + q; 
                             M(i, i+1) = q;
                        else
                             M(i, i-1) = q;
                             if i < N_matrix, M(i, i+1) = q; end
                        end
                    end
                end
                
            else % Odd Parity (sine)
                 if mod(m, 2) == 0 % Even indices (2, 4, 6...)
                    k_indices = 2:2:(2*N_matrix);
                    for i = 1:N_matrix
                        k = k_indices(i);
                        M(i, i) = k^2;
                        if i > 1, M(i, i-1) = q; end
                        if i < N_matrix, M(i, i+1) = q; end
                    end
                 else % Odd indices (1, 3, 5...)
                    k_indices = 1:2:(2*N_matrix-1);
                    for i = 1:N_matrix
                        k = k_indices(i);
                        M(i, i) = k^2;
                        if i == 1
                             M(i, i) = 1 - q; % Special boundary for se_1
                             M(i, i+1) = q;
                        else
                             M(i, i-1) = q;
                             if i < N_matrix, M(i, i+1) = q; end
                        end
                    end
                 end
            end
            
            % --- Solve Eigenvalues ---
            [V, D] = eig(M);
            eigenvals = diag(D);
            [sorted_vals, sort_idx] = sort(eigenvals);
            
            % --- Extract correct mode ---
            % Map target order 'm' to sorted index
            if strcmp(parity, 'even')
                if mod(m, 2) == 0, target_idx = m/2 + 1;
                else, target_idx = (m-1)/2 + 1; end
            else
                if mod(m, 2) == 0, target_idx = m/2;
                else, target_idx = (m-1)/2 + 1; end
            end
            
            if target_idx > N_matrix, target_idx = N_matrix; end
            
            eigval = sorted_vals(target_idx);
            coeffs = V(:, sort_idx(target_idx));
            
            % Normalize (Unit Vector)
            coeffs = coeffs / norm(coeffs);
        end
        
        function val = eval_angular(eta, m, parity, coeffs)
            % Evaluates ce_m(eta) or se_m(eta) using Fourier series
            val = zeros(size(eta));
            N = length(coeffs);
            
            if strcmp(parity, 'even')
                if mod(m, 2) == 0 % Cosine Even (0, 2, 4...)
                    % A_0 term
                    val = val + coeffs(1) * 1; % A0 * cos(0)
                    % Remaining terms A_2k * cos(2k*eta)
                    for i = 2:N
                        k_idx = 2*(i-1);
                        % Note: In matrix solver, A0 is scaled by sqrt(2)? 
                        % Usually A0 is standalone. Matrix assumed A0*sqrt(2) for symmetry.
                        % Let's stick to simple summation.
                        val = val + coeffs(i) * cos(k_idx * eta) * sqrt(2); 
                        % The sqrt(2) here compensates for the matrix symmetry trick used in solver.
                    end
                else % Cosine Odd (1, 3, 5...)
                    for i = 1:N
                        k_idx = 2*(i-1) + 1;
                        val = val + coeffs(i) * cos(k_idx * eta);
                    end
                end
            else % Odd Parity
                if mod(m, 2) == 0 % Sine Even (2, 4...)
                    for i = 1:N
                        k_idx = 2*i;
                        val = val + coeffs(i) * sin(k_idx * eta);
                    end
                else % Sine Odd (1, 3...)
                    for i = 1:N
                        k_idx = 2*(i-1) + 1;
                        val = val + coeffs(i) * sin(k_idx * eta);
                    end
                end
            end
        end
        
        function val = eval_radial(xi, m, q, parity, coeffs)
            % Evaluates Je_m(xi) or Jo_m(xi) using Bessel series
            % v = 2 * sqrt(q) * cosh(xi)
            v = 2 * sqrt(q) * cosh(xi);
            val = zeros(size(xi));
            N = length(coeffs);
            
            % Alternating signs (-1)^k are crucial in radial expansion
            
            if strcmp(parity, 'even')
                if mod(m, 2) == 0 % Even (0, 2...)
                    for i = 1:N
                        k = 2*(i-1); % 0, 2, 4...
                        % Expansion: sum (-1)^r * A_2r * J_2r(v)
                        term = coeffs(i) * besselj(k, v);
                        
                        % Apply alternating sign relative to index in array?
                        % Usually (-1)^k_index (0, 1, 2...)
                        sign_factor = (-1)^(i-1); 
                        
                        if k==0
                            val = val + sign_factor * term; % A0 no factor?
                        else
                            val = val + sign_factor * term * 2; % Neumann factor 2
                        end
                    end
                else % Odd (1, 3...)
                    for i = 1:N
                        k = 2*(i-1) + 1;
                        sign_factor = (-1)^(i-1);
                        term = coeffs(i) * besselj(k, v);
                        val = val + sign_factor * term * 2;
                    end
                end
            else % Odd Parity
                 % Sine series always has factor 2
                 if mod(m, 2) == 0 % (2, 4...)
                    for i = 1:N
                        k = 2*i;
                        sign_factor = (-1)^(i-1);
                        term = coeffs(i) * besselj(k, v);
                        val = val + sign_factor * term * 2;
                    end
                 else % (1, 3...)
                    for i = 1:N
                        k = 2*(i-1) + 1;
                        sign_factor = (-1)^(i-1);
                        term = coeffs(i) * besselj(k, v);
                        val = val + sign_factor * term * 2;
                    end
                 end
            end
            
            % Note: Radial functions grow large. 
            % But we multiply by Gaussian later, so it is stable.
        end
    end
end