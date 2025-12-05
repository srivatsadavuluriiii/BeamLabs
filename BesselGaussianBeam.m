classdef BesselGaussianBeam
    % BesselGaussianBeam: Models a Bessel function apodized by a Gaussian.
    % Reference: F. Gori, G. Guattari, and C. Padovani, "Bessel-Gauss beams," 
    % Optics Communications 64, 491-495 (1987).
    
    properties
        l               % Azimuthal topological charge (order of Bessel function)
        beta            % Transverse wave number (controls cone angle/ring spacing)
        wavelength      % Wavelength in meters
        w0              % Gaussian waist radius at z=0
        
        % Derived Physics Constants
        k               % Total wavenumber (2*pi/lambda)
        z0_gauss        % Rayleigh range of the Gaussian envelope
        alpha_cone      % Geometric cone angle (rad)
        z_max           % Maximum "diffraction-free" propagation distance
    end
    
    methods
        function obj = BesselGaussianBeam(l, beta, wavelength, w0)
            % Constructor
            % l: Topological charge (integer)
            % beta: Transverse wavenumber (rad/m). Higher beta = tighter rings.
            % w0: Width of the Gaussian apodization window.
            
            obj.l = l;
            obj.beta = beta;
            obj.wavelength = wavelength;
            obj.w0 = w0;
            
            obj.k = 2 * pi / wavelength;
            obj.z0_gauss = (pi * w0^2) / wavelength;
            
            % The geometric cone angle alpha = asin(beta/k)
            % Valid approximation for paraxial: alpha ~ beta/k
            obj.alpha_cone = asin(obj.beta / obj.k);
            
            % z_max is the geometric overlap distance where the conical 
            % waves walk off the Gaussian aperture.
            % z_max = w0 / tan(alpha) ~ k * w0 / beta
            obj.z_max = obj.k * obj.w0 / obj.beta;
        end
        
        function [field] = generate_beam_field(obj, r, phi, z, varargin)
            % Generates the complex field amplitude at distance z using
            % the Gori formulation.
            
            p_parser = inputParser;
            addRequired(p_parser, 'r');
            addRequired(p_parser, 'phi');
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            addParameter(p_parser, 'laser_linewidth_kHz', []);
            addParameter(p_parser, 'timing_jitter_ps', []);
            addParameter(p_parser, 'tx_aperture_radius', []);
            addParameter(p_parser, 'beam_tilt_x_rad', 0.0);
            addParameter(p_parser, 'beam_tilt_y_rad', 0.0);
            addParameter(p_parser, 'phase_noise_samples', []);
            addParameter(p_parser, 'symbol_time_s', []);
            
            parse(p_parser, r, phi, z, varargin{:});
            args = p_parser.Results;
            
            if ~isscalar(z)
                error('generate_beam_field expects scalar z.');
            end
            if ~isequal(size(r), size(phi))
                error('r and phi must have the same size.');
            end

            % --- Physics Calculation (Gori Eq 2.8 adapted) ---
            
            % Complex curvature parameter xi(z) = 1 + i * z / z_R
            xi = 1 + (1i * z / obj.z0_gauss);
            
            % 1. Gaussian Expansion & Amplitude Decay
            % The 1/xi term handles the Gouy phase and amplitude drop of the Gaussian envelope
            gauss_amp = (1.0 / xi) * exp(-1i * atan(z / obj.z0_gauss)); % Explicit phase extraction
            % Actually, 1/xi covers both amplitude scaling and Gouy phase shift mathematically
            term_gauss_envelope = (1.0 / xi) .* exp( -r.^2 ./ (obj.w0^2 * xi) );
            
            % 2. Bessel Component with Complex Argument
            % The argument of the Bessel function scales with propagation
            bessel_arg = obj.beta * r ./ xi;
            term_bessel = besselj(obj.l, bessel_arg);
            
            % 3. Longitudinal Corrections (The "Bessel-Gauss" axial phase correction)
            % This term describes the rapid axial phase oscillation and the envelope shift
            % Exp[ (beta^2 * w0^2 / 4) * (1/xi - 1) ] ... wait, standard form:
            % Usually expressed as exp( -i * beta^2 * z / (2k) / xi ) * exp(...)
            % Let's use the simplified compiled exponential:
            % Exp[ -i * z * beta^2 / (2k * xi) + beta^2 * w0^2 / (4 * xi) ] ??
            % Let's use the exact Gori form:
            % E = ... exp( -i * beta^2 * z / 2k ) * exp( beta^2 * z^2 / (2k * i * zR * xi) ) ...
            % Easier numerically stable form:
            val_bg_phase = exp( (obj.beta^2 * obj.w0^2 / 4) * ( (1./xi) - 1 ) );
            
            % Power Scaling
            % Normalization of BG beams is tricky as it changes with z (the core power drops).
            % We normalize so that P_tx represents the total energy in the Gaussian envelope.
            % Approx scale:
            power_scale = 1.0;
            if args.P_tx_watts > 0
                power_scale = sqrt(args.P_tx_watts);
            end
            
            % Beam Steering
            x = r .* cos(phi);
            y = r .* sin(phi);
            steering_phase = obj.k * (x * args.beam_tilt_x_rad + y * args.beam_tilt_y_rad);
            
            azimuthal_factor = exp(-1i * obj.l * phi);
            
            % Plane wave carrier phase
            propagation_phase_k = exp(1i * obj.k * z);

            % --- Noise Injection (Same as LG) ---
            phase_noise = 0.0;
            if ~isempty(args.phase_noise_samples)
                if ~isscalar(args.phase_noise_samples), error('Scalar phase noise required'); end
                phase_noise = double(args.phase_noise_samples);
            elseif ~isempty(args.laser_linewidth_kHz) && args.laser_linewidth_kHz > 0
                if isempty(args.symbol_time_s), error('symbol_time_s required'); end
                delta_nu = args.laser_linewidth_kHz * 1e3;
                sigma = sqrt(2 * pi * delta_nu * args.symbol_time_s);
                phase_noise = randn() * sigma;
            end
            
            timing_jitter_phase = 0.0;
            if ~isempty(args.timing_jitter_ps) && args.timing_jitter_ps > 0
                f_carrier = 3e8 / obj.wavelength;
                jitter_err = 2 * pi * args.timing_jitter_ps * 1e-12 * f_carrier;
                timing_jitter_phase = randn() * jitter_err;
            end
            
            noise_phasor = exp(1i * (phase_noise + timing_jitter_phase));
            
            % Combine
            field = power_scale * ...
                    term_gauss_envelope .* ...
                    term_bessel .* ...
                    val_bg_phase .* ...
                    azimuthal_factor .* ...
                    exp(1i * steering_phase) .* ...
                    propagation_phase_k .* ...
                    noise_phasor;
                
             % Aperture Mask
            if ~isempty(args.tx_aperture_radius)
                mask = double(r <= args.tx_aperture_radius);
                field = field .* mask;
            end
        end
        
        function intensity = calculate_intensity(obj, r, phi, z, varargin)
            field = obj.generate_beam_field(r, phi, z, varargin{:});
            intensity = abs(field).^2;
        end
        
        function pn = generate_phase_noise_sequence(~, num_symbols, symbol_time_s, laser_linewidth_kHz, seed)
            if nargin > 4 && ~isempty(seed), rng(seed); end
            if isempty(laser_linewidth_kHz) || laser_linewidth_kHz <= 0
                pn = zeros(1, num_symbols); return;
            end
            delta_nu = laser_linewidth_kHz * 1e3;
            phase_variance = 2 * pi * delta_nu * symbol_time_s;
            increments = randn(1, num_symbols) * sqrt(phase_variance);
            pn = cumsum(increments);
        end
        
        function params = get_beam_parameters(obj, z)
            % Returns parameters specific to Bessel propagation
            params.z = z;
            params.in_bessel_zone = (z < obj.z_max);
            params.fraction_zmax = z / obj.z_max;
            
            % Gaussian envelope waist at z (The Bessel core sits inside this)
            params.envelope_waist_z = obj.w0 * sqrt(1 + (z/obj.z0_gauss)^2);
            
            % Approximate core radius (central lobe size)
            % For J_l(x), first root is roughly at x ~ 2.4 (for l=0).
            % Argument is beta * r. Core radius ~ 2.405 / beta
            % This stays CONSTANT in the Bessel zone.
            params.core_radius_approx = 2.4048 / obj.beta; 
            
            params.z_max = obj.z_max;
        end
        
        function propagation_summary(obj, z_distances)
            fprintf('\nPropagation Summary for BG beam (l=%d, beta=%.1f rad/m):\n', obj.l, obj.beta);
            fprintf('z_max (Diffraction Free Range) = %.2f m\n', obj.z_max);
            fprintf('%s\n', repmat('-', 1, 85));
            fprintf('%-10s | %-12s | %-15s | %-15s | %-10s\n', ...
                'Dist (m)', 'Region', 'Envelope w(z)', 'Core Size', 'On-Axis Int');
            fprintf('%s\n', repmat('-', 1, 85));
            
            for i = 1:length(z_distances)
                z = z_distances(i);
                p = obj.get_beam_parameters(z);
                
                if p.in_bessel_zone
                    region = 'Bessel Zone';
                    % In Bessel zone, core size is constant determined by beta
                    core_str = sprintf('~%.3f mm (C)', p.core_radius_approx*1e3);
                else
                    region = 'Far Field';
                    % In far field, it becomes a ring, core dissolves
                    core_str = 'Diffracting';
                end
                
                % Calculate on-axis (or near axis for l>0) intensity for checking
                % For l=0, r=0. For l>0, check r=small epsilon
                if obj.l == 0
                    int_val = obj.calculate_intensity(0, 0, z);
                else
                    % Look slightly off axis to find peak of ring
                    r_peak = p.core_radius_approx; 
                    int_val = obj.calculate_intensity(r_peak, 0, z);
                end
                
                fprintf('%8.1f   | %-12s | %8.2f mm      | %-15s | %.2e\n', ...
                    z, region, p.envelope_waist_z*1e3, core_str, int_val);
            end
        end
        
        function summary = get_tx_parameters_summary(obj, varargin)
            % Same utility as LG beam for link budget summaries
            p_parser = inputParser;
            addParameter(p_parser, 'P_tx_watts', 1.0);
            addParameter(p_parser, 'laser_linewidth_kHz', []);
            addParameter(p_parser, 'timing_jitter_ps', []);
            addParameter(p_parser, 'tx_aperture_radius', []);
            addParameter(p_parser, 'beam_tilt_x_rad', 0.0);
            addParameter(p_parser, 'beam_tilt_y_rad', 0.0);
            parse(p_parser, varargin{:});
            args = p_parser.Results;
            
            summary.P_tx_watts = args.P_tx_watts;
            summary.P_tx_dBm = 10 * log10(args.P_tx_watts * 1000);
            summary.beta_rad_m = obj.beta;
            summary.cone_angle_deg = rad2deg(obj.alpha_cone);
            summary.diffraction_free_range_m = obj.z_max;
            
            % ... (Jitter/Phase noise calcs same as LG) ...
             if ~isempty(args.laser_linewidth_kHz)
                summary.laser_linewidth_kHz = args.laser_linewidth_kHz;
            else
                summary.laser_linewidth_kHz = 0.0;
            end
            
            if ~isempty(args.timing_jitter_ps)
                summary.timing_jitter_ps = args.timing_jitter_ps;
            else
                summary.timing_jitter_ps = 0.0;
            end
            
            if summary.timing_jitter_ps > 0
                f_carrier = 3e8 / obj.wavelength;
                rms_rad = 2 * pi * summary.timing_jitter_ps * 1e-12 * f_carrier;
                summary.timing_jitter_phase_rms_rad = rms_rad;
                summary.timing_jitter_phase_rms_deg = rad2deg(rms_rad);
            else
                summary.timing_jitter_phase_rms_rad = 0;
                summary.timing_jitter_phase_rms_deg = 0;
            end
        end
    end
end