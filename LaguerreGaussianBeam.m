classdef LaguerreGaussianBeam
    % LaguerreGaussianBeam: A class to model and analyze LG beams.
    % Includes support for phase noise, timing jitter, and aperture loss.
    
    properties
        p               % Radial mode index (Integer)
        l               % Azimuthal mode index (Integer)
        wavelength      % Wavelength in meters
        w0              % Beam waist radius at z=0 in meters
        k               % Wavenumber
        M_squared       % Beam quality factor
        z0_fundamental  % Rayleigh range (fundamental mode)
        z_R             % Rayleigh range (actual mode)
        C_norm          % Normalization constant
    end
    
    methods
        function obj = LaguerreGaussianBeam(p, l, wavelength, w0)
            if p < 0
                error('Radial index p must be non-negative, got p=%d', p);
            end
            
            obj.p = p;
            obj.l = l;
            obj.wavelength = wavelength;
            obj.w0 = w0;
            
            % Derived constants
            obj.k = 2 * pi / wavelength;
            obj.M_squared = 2 * p + abs(l) + 1;
            obj.z0_fundamental = (pi * w0^2) / wavelength;
            obj.z_R = obj.z0_fundamental / obj.M_squared;
            
            % Normalization constant
            % Note: using exp(gammaln) handles large factorials better
            fact_p = exp(gammaln(p + 1));
            fact_pl = exp(gammaln(p + abs(l) + 1));
            obj.C_norm = sqrt(2.0 * fact_p / (pi * fact_pl));
        end
        
        function val = beam_waist(obj, z)
            val = obj.w0 * sqrt(1 + ( (obj.M_squared * z) / obj.z0_fundamental ).^2);
        end
        
        function val = radius_of_curvature(obj, z)
            % Handle scalar or array inputs
            val = zeros(size(z));
            mask = abs(z) >= 1e-12;
            
            % R(z) = z * (1 + (z_R/z)^2)
            val(mask) = z(mask) .* (1 + (obj.z_R ./ z(mask)).^2);
            
            % Set near-zero z to Infinity
            val(~mask) = inf;
        end
        
        function val = gouy_phase(obj, z)
            val = obj.M_squared * atan(z / obj.z0_fundamental);
        end
        
        function [theta_0, theta_eff] = effective_divergence_angle(obj)
            theta_0 = obj.wavelength / (pi * obj.w0);
            theta_eff = obj.M_squared * theta_0;
        end
        
        function field = generate_beam_field(obj, r, phi, z, varargin)
            % FIX: Use a local variable 'p_parser', DO NOT use 'obj.p'
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
            
            % Ensure r and phi are same size (broadcasting check)
            if ~isequal(size(r), size(phi))
                error('r and phi must have the same size.');
            end
            
            w_z = obj.beam_waist(z);
            R_z = obj.radius_of_curvature(z);
            psi_z = obj.gouy_phase(z);
            
            arg_L = 2 * r.^2 / w_z^2;
            
            % Numerical Generalized Laguerre Polynomial
            % obj.p is now correctly the integer index, not a parser
            L_p_l = LaguerreGaussianBeam.numerical_genlaguerre(obj.p, abs(obj.l), arg_L);
            
            power_scale = 0;
            if args.P_tx_watts > 0
                power_scale = sqrt(args.P_tx_watts);
            end
            
            amplitude_factor = obj.C_norm * (1.0 / w_z) * power_scale;
            
            radial_factor = ( (sqrt(2)*r./w_z).^abs(obj.l) ) .* ...
                            L_p_l .* exp(-(r.^2) / w_z^2);
            
            % Aperture Mask
            mask_factor = 1.0;
            if ~isempty(args.tx_aperture_radius)
                mask_factor = double(r <= args.tx_aperture_radius);
            end
            
            % Beam Steering
            x = r .* cos(phi);
            y = r .* sin(phi);
            steering_phase = obj.k * (x * args.beam_tilt_x_rad + y * args.beam_tilt_y_rad);
            
            azimuthal_factor = exp(-1i * obj.l * phi);
            
            if isinf(R_z)
                curvature_factor = 1.0;
            else
                curvature_phase = -1i * obj.k * r.^2 / (2 * R_z);
                curvature_factor = exp(curvature_phase);
            end
            
            gouy_term = exp(-1i * psi_z);
            
            % Noise logic
            phase_noise = 0.0;
            if ~isempty(args.phase_noise_samples)
                if ~isscalar(args.phase_noise_samples)
                   error('phase_noise_samples must be scalar.');
                end
                phase_noise = double(args.phase_noise_samples);
            elseif ~isempty(args.laser_linewidth_kHz) && args.laser_linewidth_kHz > 0
                if isempty(args.symbol_time_s)
                    error('symbol_time_s required for on-the-fly phase noise.');
                end
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
            
            propagation_phase = exp(1i * (obj.k * z + phase_noise + timing_jitter_phase));
            
            % Combine terms
            field = amplitude_factor * radial_factor .* azimuthal_factor .* ...
                    curvature_factor .* exp(1i * steering_phase) * ...
                    gouy_term * propagation_phase;
            
            if ~isempty(args.tx_aperture_radius)
                field = field .* mask_factor;
            end
        end
        
        function intensity = calculate_intensity(obj, r, phi, z, varargin)
             field = obj.generate_beam_field(r, phi, z, varargin{:});
             intensity = abs(field).^2;
        end
        
        function pn = generate_phase_noise_sequence(~, num_symbols, symbol_time_s, laser_linewidth_kHz, seed)
            if nargin > 4 && ~isempty(seed)
                rng(seed);
            end
            
            if isempty(laser_linewidth_kHz) || laser_linewidth_kHz <= 0
                pn = zeros(1, num_symbols);
                return;
            end
            
            delta_nu = laser_linewidth_kHz * 1e3;
            phase_variance = 2 * pi * delta_nu * symbol_time_s;
            
            increments = randn(1, num_symbols) * sqrt(phase_variance);
            pn = cumsum(increments);
        end
        
        function transmission = calculate_tx_aperture_loss(obj, tx_aperture_radius, varargin)
            % Defaults
            r_max_factor = 8.0;
            n_r = 5000;
            if nargin > 2
                % Simple parsing for positional optional args if strictly needed
                % keeping it simple for now
            end
            
            if isempty(tx_aperture_radius) || tx_aperture_radius <= 0
                transmission = 1.0;
                return;
            end
            
            r_max = max(tx_aperture_radius * r_max_factor, obj.w0 * r_max_factor);
            r_array = linspace(0, r_max, n_r);
            
            % Intensity at z=0, phi=0
            phi_dummy = zeros(size(r_array));
            intensity = obj.calculate_intensity(r_array, phi_dummy, 0);
            
            integrand = intensity .* 2 * pi .* r_array;
            total_power = trapz(r_array, integrand);
            
            mask = r_array <= tx_aperture_radius;
            r_ap = r_array(mask);
            int_ap = integrand(mask);
            
            if isempty(r_ap)
                transmission = 0.0;
                return;
            end
            
            aperture_power = trapz(r_ap, int_ap);
            
            if total_power > 0
                transmission = aperture_power / total_power;
            else
                transmission = 0.0;
            end
        end
        
        function overlap = overlap_with(obj, other_beam, r_max_factor, n_r, n_phi)
            if nargin < 3, r_max_factor = 6.0; end
            if nargin < 4, n_r = 800; end
            if nargin < 5, n_phi = 360; end
            
            r_max = max(obj.w0, other_beam.w0) * r_max_factor;
            
            % 1D arrays
            r_vec = linspace(0, r_max, n_r);
            phi_vec = linspace(0, 2*pi, n_phi + 1); 
            phi_vec(end) = []; 
            
            [R, PHI] = meshgrid(r_vec, phi_vec);
            
            % Generate fields
            E1 = obj.generate_beam_field(R, PHI, 0.0, 'P_tx_watts', 1.0);
            E2 = other_beam.generate_beam_field(R, PHI, 0.0, 'P_tx_watts', 1.0);
            
            integrand = conj(E1) .* E2 .* R;
            
            int_over_phi = trapz(phi_vec, integrand, 1);
            overlap = trapz(r_vec, int_over_phi);
        end
        
        function params = get_beam_parameters(obj, z)
            params.z = z;
            params.w_z = obj.beam_waist(z);
            params.R_z = obj.radius_of_curvature(z);
            params.gouy_phase = obj.gouy_phase(z);
            params.beam_expansion = params.w_z / obj.w0;
        end
        
        function propagation_summary(obj, z_distances)
            fprintf('\nPropagation Summary for LG_%d^%d beam:\n', obj.p, obj.l);
            fprintf('Distance (m) | w(z) (mm) | R(z) (m) | Gouy Phase (rad) | w(z)/w0\n');
            fprintf('%s\n', repmat('-', 1, 65));
            
            for i = 1:length(z_distances)
                z = z_distances(i);
                prms = obj.get_beam_parameters(z);
                
                if isinf(prms.R_z)
                    R_str = 'Inf';
                else
                    R_str = sprintf('%.2f', prms.R_z);
                end
                
                fprintf('%8.1f   | %7.2f   | %8s | %11.3f    | %.2f\n', ...
                    z, prms.w_z*1e3, R_str, prms.gouy_phase, prms.beam_expansion);
            end
        end
        
        function summary = get_tx_parameters_summary(obj, varargin)
            % FIX: Use a local variable 'p_parser', DO NOT use 'obj.p'
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
            if args.P_tx_watts > 0
                summary.P_tx_dBm = 10 * log10(args.P_tx_watts * 1000);
            else
                summary.P_tx_dBm = -inf;
            end
            
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
            
            summary.tx_aperture_radius_m = args.tx_aperture_radius;
            summary.beam_tilt_x_rad = args.beam_tilt_x_rad;
            summary.beam_tilt_y_rad = args.beam_tilt_y_rad;
            summary.beam_tilt_x_deg = rad2deg(args.beam_tilt_x_rad);
            summary.beam_tilt_y_deg = rad2deg(args.beam_tilt_y_rad);
            
            if ~isempty(args.tx_aperture_radius)
                trans = obj.calculate_tx_aperture_loss(args.tx_aperture_radius);
                summary.tx_aperture_transmission = trans;
                if trans > 0
                    summary.tx_aperture_loss_dB = -10 * log10(trans);
                else
                    summary.tx_aperture_loss_dB = inf;
                end
            else
                summary.tx_aperture_transmission = 1.0;
                summary.tx_aperture_loss_dB = 0.0;
            end
            
            % Jitter calculations
            if summary.timing_jitter_ps > 0
                f_carrier = 3e8 / obj.wavelength;
                rms_rad = 2 * pi * summary.timing_jitter_ps * 1e-12 * f_carrier;
                summary.timing_jitter_phase_rms_rad = rms_rad;
                summary.timing_jitter_phase_rms_deg = rad2deg(rms_rad);
                summary.timing_jitter_phase_rms_cycles = rms_rad / (2*pi);
            else
                summary.timing_jitter_phase_rms_rad = 0;
                summary.timing_jitter_phase_rms_deg = 0;
                summary.timing_jitter_phase_rms_cycles = 0;
            end
            
            % Phase noise calculations (illustrative)
            if summary.laser_linewidth_kHz > 0
                symbol_rate = 1e9;
                sym_time = 1.0/symbol_rate;
                pn_rms = sqrt(2*pi*summary.laser_linewidth_kHz*1e3*sym_time);
                summary.phase_noise_rms_rad = pn_rms;
                summary.phase_noise_rms_deg = rad2deg(pn_rms);
            else
                summary.phase_noise_rms_rad = 0;
                summary.phase_noise_rms_deg = 0;
            end
        end
    end
    
    methods (Static)
        function val = numerical_genlaguerre(n, alpha, x)
            % Evaluates Generalized Laguerre Polynomial L_n^alpha(x)
            
            if n == 0
                val = ones(size(x));
                return;
            end
            if n == 1
                val = 1 + alpha - x;
                return;
            end
            
            L_k_minus_1 = ones(size(x));      % L_0
            L_k = 1 + alpha - x;              % L_1
            
            for k = 1:(n-1)
                % Recurrence
                % L_next = [ (2k + 1 + alpha - x)L_k - (k + alpha)L_prev ] / (k+1)
                L_next = ( (2*k + 1 + alpha - x) .* L_k - (k + alpha) * L_k_minus_1 ) / (k + 1);
                
                L_k_minus_1 = L_k;
                L_k = L_next;
            end
            val = L_k;
        end
    end
end