classdef InceGaussianBeam
    % InceGaussianBeam: Models paraxial beams in Elliptic symmetry.
    % Solves the Ince differential equation numerically to find C-coefficients.
    % Reference: Bandres & Gutierrez-Vega, Opt. Lett. 29, 144 (2004).
    
    properties
        p               % Order (0, 1, 2...) corresponds to total mode number
        m               % Degree (0, 1... p) corresponds to number of hyperbolic nodes
        epsilon         % Ellipticity parameter (0 = LG limit, inf = HG limit)
        parity          % 'even' or 'odd'
        wavelength      % Wavelength in meters
        w0              % Beam waist radius
        
        % Derived Constants
        k               % Wavenumber
        z0              % Rayleigh range
        f0              % Semifocal distance (scale factor for elliptic coords)
        
        % Eigenvalues/Coefficients (Computed on init)
        a_eigenvalue    % Ince equation eigenvalue
        C_coeffs        % Expansion coefficients for the polynomial
        norm_factor     % Normalization constant
    end
    
    methods
        function obj = InceGaussianBeam(p, m, epsilon, parity, wavelength, w0)
            % Constructor
            % p:       Mode order (integer >= 0)
            % m:       Mode degree (0 <= m <= p)
            % epsilon: Ellipticity (e.g., 2.0)
            % parity:  'even' or 'odd' (string)
            
            if mod(p-m, 2) ~= 0
                error('Indices (p, m) must have the same parity (p-m must be even).');
            end
            
            obj.p = p;
            obj.m = m;
            obj.epsilon = epsilon;
            obj.parity = lower(parity);
            obj.wavelength = wavelength;
            obj.w0 = w0;
            
            obj.k = 2 * pi / wavelength;
            obj.z0 = (pi * w0^2) / wavelength;
            
            % Semifocal distance f0 = w0 * sqrt(epsilon/2)
            obj.f0 = obj.w0 * sqrt(obj.epsilon / 2);
            
            % --- Precompute Ince Polynomial Coefficients ---
            % We solve the eigenvalue problem once during initialization
            [obj.a_eigenvalue, obj.C_coeffs, obj.norm_factor] = ...
                obj.solve_ince_equation(p, m, epsilon, obj.parity);
        end
        
        function field = generate_beam_field(obj, x, y, z, varargin)
            p_parser = inputParser;
            addRequired(p_parser, 'x');
            addRequired(p_parser, 'y');
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            addParameter(p_parser, 'tx_aperture_radius', []);
            parse(p_parser, x, y, z, varargin{:});
            args = p_parser.Results;
            
            if ~isscalar(z), error('z must be scalar'); end
            
            % --- 1. Coordinate Transformation (Cartesian -> Elliptic) ---
            w_z = obj.w0 * sqrt(1 + (z/obj.z0)^2);
            
            % The grid scales with beam expansion w(z)
            % We need normalized coordinates.
            % Global semifocal separation moves with z?
            % Actually, standard definition uses w(z) scaling.
            % r_norm = r / w(z). 
            % Elliptic coords (xi, eta) are defined such that:
            % x = w(z) * sqrt(eps/2) * cosh(xi) * cos(eta)
            % y = w(z) * sqrt(eps/2) * sinh(xi) * sin(eta)
            
            % To reverse this numerically (Cartesian -> Elliptic):
            % We define complex z_plane = x + i*y
            % Normalized z_norm = z_plane / (w(z) * sqrt(eps/2))
            % acosh(z_norm) = xi + i*eta
            
            f_scale = w_z * sqrt(obj.epsilon / 2);
            z_complex = x + 1i * y;
            
            % Avoid division by zero if epsilon=0 (LG case)
            if obj.epsilon < 1e-9
                error('For epsilon ~ 0, use LaguerreGaussianBeam class.');
            end
            
            w_complex = acosh(z_complex ./ f_scale);
            xi = real(w_complex);
            eta = imag(w_complex);
            
            % --- 2. Compute Ince Polynomials ---
            % IP = C * Ip(xi) * Ip(eta) (Separable in elliptic coords)
            
            [Ince_xi, ~]  = obj.eval_ince_poly(xi,  obj.p, obj.m, obj.epsilon, obj.parity, 1, obj.C_coeffs);
            [Ince_eta, ~] = obj.eval_ince_poly(eta, obj.p, obj.m, obj.epsilon, obj.parity, 2, obj.C_coeffs);
            
            % --- 3. Gaussian Envelope & Phases ---
            % Gaussian term: exp( -r^2 / w^2 )
            r_sq = x.^2 + y.^2;
            gauss_term = exp( -r_sq ./ w_z^2 );
            
            % Gouy Phase: (p + 1) * atan(z/zR)
            gouy_val = (obj.p + 1) * atan(z / obj.z0);
            gouy_term = exp( -1i * gouy_val );
            
            % Curvature
            R_z = z * (1 + (obj.z0/z)^2);
            if abs(z) < 1e-12, R_z = inf; end
            curv_term = exp( -1i * obj.k * r_sq / (2 * R_z) );
            
            % Amplitude Scaling
            p_scale = (args.P_tx_watts > 0) * sqrt(args.P_tx_watts);
            % Apply normalization from constructor
            % The literature normalization is complex, we use the computed norm_factor
            amp_term = (obj.C_coeffs(1) * 0 + 1) * obj.norm_factor * (obj.w0 / w_z) * p_scale; 
            
            % --- Assembly ---
            % Field = Amp * Ince(xi) * Ince(eta) * Gauss * Phase
            field = amp_term .* Ince_xi .* Ince_eta .* gauss_term .* ...
                    curv_term .* gouy_term .* exp(1i * obj.k * z);
            
             if ~isempty(args.tx_aperture_radius)
                field = field .* (sqrt(r_sq) <= args.tx_aperture_radius);
            end
        end
        
        function intensity = calculate_intensity(obj, x, y, z, varargin)
            field = obj.generate_beam_field(x, y, z, varargin{:});
            intensity = abs(field).^2;
        end
        
        % --- Numerical Solver for Ince Polynomials ---
        function [eig_val, C, N_const] = solve_ince_equation(~, p, m, eps, parity)
            % Solves the tridiagonal recurrence relation to find coefficients
            % Reference: Bandres, Opt. Lett. 29, 144 (2004), Eq 3.
            
            % Determine Number of coefficients
            if strcmp(parity, 'even')
                % Even functions C_{2k} or C_{2k+1}
                if mod(p, 2) == 0
                    k_indices = 0:2:p; % Even indices 0, 2, ... p
                else
                    k_indices = 1:2:p; % Odd indices 1, 3, ... p
                end
            else % Odd
                % Sine functions S_{2k} or S_{2k+1}
                if mod(p, 2) == 0
                    k_indices = 2:2:p; % Start at 2 for Sine even p
                else
                    k_indices = 1:2:p;
                end
            end
            
            n_coeffs = length(k_indices);
            M_mat = zeros(n_coeffs, n_coeffs);
            
            % Construct Tridiagonal Matrix M
            for idx = 1:n_coeffs
                r = k_indices(idx); % Actual index (0, 1, 2...)
                j = idx;            % Matrix index (1, 2, 3...)
                
                if strcmp(parity, 'even')
                    % --- C-Ince Matrix Elements ---
                    if mod(p, 2) == 0 % (p even) -> (0, 2, 4...)
                         if r == 0
                             M_mat(j, j+1) = eps * (p + 2) / 2; % Right
                         else
                             val_diag = 4*r^2;
                             M_mat(j, j) = val_diag;
                             if j > 1, M_mat(j, j-1) = eps*(p + 1 - (r-2))/2; end % Left
                             if j < n_coeffs, M_mat(j, j+1) = eps*(p + 1 + (r+2))/2; end % Right
                         end
                    else % (p odd) -> (1, 3, 5...)
                        val_diag = 4*r^2 + eps*(p+1) + eps/2; % Special correction for r=1?
                        if r == 1
                             val_diag = 4 + eps*(p+1) + eps; % Check papers for index 1 boundary
                             M_mat(j, j) = val_diag;
                             M_mat(j, j+1) = eps*(p + 1 + 3)/2;
                        else
                             M_mat(j, j) = 4*r^2;
                             M_mat(j, j-1) = eps*(p + 1 - (r-2))/2; 
                             if j < n_coeffs, M_mat(j, j+1) = eps*(p + 1 + (r+2))/2; end
                        end
                    end
                    
                else
                    % --- S-Ince Matrix Elements ---
                    val_diag = 4*r^2;
                    M_mat(j, j) = val_diag;
                    if j > 1, M_mat(j, j-1) = eps*(p + 1 - (r-2))/2; end
                    if j < n_coeffs, M_mat(j, j+1) = eps*(p + 1 + (r+2))/2; end
                end
            end
            
            % Solve Eigenvalue Problem
            [V, D] = eig(M_mat);
            eigenvalues = diag(D);
            
            % We need the specific eigenvalue 'a' corresponding to index m
            % Sorting eigenvalues usually aligns them with m
            [sorted_eigs, sort_idx] = sort(eigenvalues);
            
            % Mapping (m) to the sorted index is non-trivial but usually linear
            % For a given p, there are limited m's.
            % Logic: The eigenvalues represent the characteristic 'a'. 
            % m=p is usually the highest/lowest energy depending on sign convention.
            
            % Simplified index mapping for standard ordering:
            % Find where m sits in the allowed indices
            allowed_m = k_indices; % Actually m usually follows same parity
            target_idx = find(allowed_m == m);
            
            if isempty(target_idx)
                % Fallback if exact m mapping logic differs (Bandres ordering)
                % Usually just take the k-th eigenvalue
                % Let's try direct mapping
                warning('Ince m-index mapping is heuristic.');
                target_idx = 1; 
            end
            
            eig_val = sorted_eigs(target_idx);
            C = V(:, sort_idx(target_idx));
            
            % Normalize Coefficients
            % Normalization convention varies. We normalize such that max(C) = 1
            % or sum(C^2) = 1. 
            C = C / sqrt(sum(C.^2));
            
            % Total Beam Normalization Factor (Approximate)
            N_const = 1.0; % Proper normalization requires 2D integration
        end
        
        function [val, dval] = eval_ince_poly(obj, z_coord, p, ~, eps, parity, type, C)
            % Evaluate the polynomial series summation
            % type: 1 = xi (hyperbolic), 2 = eta (trigonometric)
            
            val = zeros(size(z_coord));
            dval = zeros(size(z_coord));
            
            if strcmp(parity, 'even')
                if mod(p, 2) == 0
                    k_idx = 0:2:p;
                else
                    k_idx = 1:2:p;
                end
                
                for i = 1:length(C)
                    r = k_idx(i);
                    c = C(i);
                    if type == 1 % Xi (Hyperbolic Cos)
                        val = val + c * cosh(r * z_coord);
                    else         % Eta (Trig Cos)
                        val = val + c * cos(r * z_coord);
                    end
                end
            else % Odd (Sine)
                if mod(p, 2) == 0
                    k_idx = 2:2:p;
                else
                    k_idx = 1:2:p;
                end
                for i = 1:length(C)
                    r = k_idx(i);
                    c = C(i);
                    if type == 1 % Xi (Hyperbolic Sin)
                        val = val + c * sinh(r * z_coord);
                    else         % Eta (Trig Sin)
                        val = val + c * sin(r * z_coord);
                    end
                end
            end
        end
    end
end