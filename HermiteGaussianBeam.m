classdef HermiteGaussianBeam
    % HermiteGaussianBeam: Models paraxial beams in Cartesian symmetry.
    % Defined by indices (n, m) corresponding to x and y axes.
    
    properties
        n               % Mode index along x-axis (integer >= 0)
        m               % Mode index along y-axis (integer >= 0)
        wavelength      % Wavelength in meters
        w0              % Beam waist radius at z=0
        
        % Derived Constants
        k               % Wavenumber
        z0              % Rayleigh range
        C_norm          % Normalization constant
    end
    
    methods
        function obj = HermiteGaussianBeam(n, m, wavelength, w0)
            if n < 0 || m < 0
                error('Indices n and m must be non-negative.');
            end
            
            obj.n = n;
            obj.m = m;
            obj.wavelength = wavelength;
            obj.w0 = w0;
            
            obj.k = 2 * pi / wavelength;
            obj.z0 = (pi * w0^2) / wavelength;
            
            % Normalization Factor for HG_nm
            % C_nm = sqrt( 2 / (pi * w0^2 * 2^(n+m) * n! * m!) )
            % We use gammaln to handle large factorials safely
            log_denom = log(pi) + 2*log(w0) + (n+m)*log(2) + gammaln(n+1) + gammaln(m+1);
            obj.C_norm = sqrt(2 * exp(-log_denom)); 
            % Note: The '2' in numerator and w0^2 in denom depends on 
            % whether you normalize over integral |E|^2 dx dy = 1.
            % This factor ensures unitary power.
        end
        
        function val = beam_waist(obj, z)
            val = obj.w0 * sqrt(1 + (z / obj.z0).^2);
        end
        
        function val = radius_of_curvature(obj, z)
            val = zeros(size(z));
            mask = abs(z) >= 1e-12;
            val(mask) = z(mask) .* (1 + (obj.z0 ./ z(mask)).^2);
            val(~mask) = inf;
        end
        
        function val = gouy_phase(obj, z)
            % HG Gouy phase depends on sum of indices (N = n + m)
            % Phi(z) = (n + m + 1) * atan(z/z0)
            val = (obj.n + obj.m + 1) * atan(z / obj.z0);
        end
        
        function field = generate_beam_field(obj, x_in, y_in, z, varargin)
            % Inputs can be (x, y) Cartesian OR (r, phi) Cylindrical.
            % To maintain polymorphism with LG/BG codes, we detect input type.
            % If inputs are named 'r' and 'phi' in the calling script, 
            % we assume they are polar and convert.
            
            p_parser = inputParser;
            addRequired(p_parser, 'in1'); % x or r
            addRequired(p_parser, 'in2'); % y or phi
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            addParameter(p_parser, 'laser_linewidth_kHz', []);
            addParameter(p_parser, 'timing_jitter_ps', []);
            addParameter(p_parser, 'tx_aperture_radius', []);
            addParameter(p_parser, 'beam_tilt_x_rad', 0.0);
            addParameter(p_parser, 'beam_tilt_y_rad', 0.0);
            addParameter(p_parser, 'phase_noise_samples', []);
            addParameter(p_parser, 'symbol_time_s', []);
            addParameter(p_parser, 'is_polar_input', false); % specific flag
            
            parse(p_parser, x_in, y_in, z, varargin{:});
            args = p_parser.Results;
            
            % Coordinate conversion logic
            if args.is_polar_input
                % Treat in1 as r, in2 as phi
                r = x_in; phi = y_in;
                x = r .* cos(phi);
                y = r .* sin(phi);
            else
                % Treat as Cartesian x, y
                x = x_in;
                y = y_in;
                % Calc r for aperture/curvature logic
                r = sqrt(x.^2 + y.^2); 
            end
            
            if ~isscalar(z), error('z must be scalar'); end
            
            % --- Beam Parameters ---
            w_z = obj.beam_waist(z);
            R_z = obj.radius_of_curvature(z);
            psi_z = obj.gouy_phase(z);
            
            % --- Hermite Polynomials ---
            % Arguments are sqrt(2)*x / w(z)
            scale_factor = sqrt(2) / w_z;
            H_n = HermiteGaussianBeam.numerical_hermite(obj.n, scale_factor * x);
            H_m = HermiteGaussianBeam.numerical_hermite(obj.m, scale_factor * y);
            
            % --- Amplitude Term ---
            % E0 * (w0/wz) * ...
            power_scale = (args.P_tx_watts > 0) * sqrt(args.P_tx_watts);
            amp_term = obj.C_norm * (obj.w0 / w_z) * power_scale;
            
            % --- Gaussian Envelope ---
            gauss_env = exp( - (x.^2 + y.^2) / w_z^2 );
            
            % --- Phase Terms ---
            if isinf(R_z)
                curv_phase = 1.0;
            else
                curv_phase = exp( -1i * obj.k * (x.^2 + y.^2) / (2 * R_z) );
            end
            
            gouy_term = exp(1i * psi_z); % Definition convention vary, usually exp(i*psi)
            
            % --- Steering & Noise ---
            steering = exp(1i * obj.k * (x * args.beam_tilt_x_rad + y * args.beam_tilt_y_rad));
            
            phase_noise = 0;
            if ~isempty(args.laser_linewidth_kHz) && args.laser_linewidth_kHz > 0
                 delta_nu = args.laser_linewidth_kHz * 1e3;
                 sigma = sqrt(2 * pi * delta_nu * args.symbol_time_s);
                 phase_noise = randn() * sigma;
            end
            
            jitter_phase = 0;
            if ~isempty(args.timing_jitter_ps) && args.timing_jitter_ps > 0
                f_c = 3e8 / obj.wavelength;
                rms = 2*pi*args.timing_jitter_ps*1e-12*f_c;
                jitter_phase = randn() * rms;
            end
            
            prop_phase = exp(1i * (obj.k * z + phase_noise + jitter_phase));
            
            % --- Assembly ---
            field = amp_term .* ...
                    H_n .* H_m .* ...
                    gauss_env .* ...
                    curv_phase .* ...
                    gouy_term .* ...
                    steering .* ...
                    prop_phase;
            
            if ~isempty(args.tx_aperture_radius)
                mask = (r <= args.tx_aperture_radius);
                field = field .* mask;
            end
        end
        
        function intensity = calculate_intensity(obj, in1, in2, z, varargin)
            field = obj.generate_beam_field(in1, in2, z, varargin{:});
            intensity = abs(field).^2;
        end
        
        function params = get_beam_parameters(obj, z)
            params.z = z;
            params.w_z = obj.beam_waist(z);
            params.R_z = obj.radius_of_curvature(z);
            params.gouy_phase = obj.gouy_phase(z);
            params.M_squared_x = 2*obj.n + 1;
            params.M_squared_y = 2*obj.m + 1;
        end
        
        function summary = get_tx_parameters_summary(obj, varargin)
             % Re-use the robust logic from previous classes
             p_parser = inputParser;
             addParameter(p_parser, 'P_tx_watts', 1.0);
             % ... add other params as needed ...
             parse(p_parser, varargin{:});
             summary = p_parser.Results;
             summary.type = sprintf('HG %d %d', obj.n, obj.m);
        end
    end
    
    methods (Static)
        function H_val = numerical_hermite(n, x)
            % Evaluates Physicists' Hermite Polynomial H_n(x)
            % H_0 = 1
            % H_1 = 2x
            % H_{n+1} = 2x H_n - 2n H_{n-1}
            
            if n == 0
                H_val = ones(size(x));
                return;
            end
            if n == 1
                H_val = 2 * x;
                return;
            end
            
            H_km1 = ones(size(x));  % H_0
            H_k = 2 * x;            % H_1
            
            for k = 1:(n-1)
                % Recurrence relation
                H_next = 2 * x .* H_k - 2 * k * H_km1;
                
                H_km1 = H_k;
                H_k = H_next;
            end
            H_val = H_k;
        end
    end
end