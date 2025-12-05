classdef GaussianBeam
    % GaussianBeam: Models the fundamental TEM00 paraxial Gaussian beam.
    % This is the baseline for optical communications and laser physics.
    
    properties
        wavelength      % Wavelength in meters
        w0              % Beam waist radius at z=0 (1/e^2 intensity)
        
        % Derived Constants
        k               % Wavenumber
        z0              % Rayleigh range (diffraction length)
        theta_div       % Far-field divergence angle (half-angle)
    end
    
    methods
        function obj = GaussianBeam(wavelength, w0)
            % Constructor
            obj.wavelength = wavelength;
            obj.w0 = w0;
            
            % Fundamental constants
            obj.k = 2 * pi / wavelength;
            obj.z0 = (pi * w0^2) / wavelength;
            obj.theta_div = wavelength / (pi * w0);
        end
        
        function val = beam_waist(obj, z)
            % w(z) = w0 * sqrt(1 + (z/z0)^2)
            val = obj.w0 * sqrt(1 + (z / obj.z0).^2);
        end
        
        function val = radius_of_curvature(obj, z)
            % R(z) = z * (1 + (z0/z)^2)
            val = zeros(size(z));
            mask = abs(z) >= 1e-12;
            val(mask) = z(mask) .* (1 + (obj.z0 ./ z(mask)).^2);
            val(~mask) = inf; % Planar wavefront at waist
        end
        
        function val = gouy_phase(obj, z)
            % Psi(z) = atan(z/z0)
            % Note: Fundamental mode has factor of 1, not (N+1)
            val = atan(z / obj.z0);
        end
        
        function field = generate_beam_field(obj, r, phi, z, varargin)
            % Generates the complex scalar field E(r,phi,z).
            % Accepts polar inputs (r, phi) by default to match LG/BG interfaces.
            
            p_parser = inputParser;
            addRequired(p_parser, 'r');
            addRequired(p_parser, 'phi'); % Kept for interface consistency, though unused for TEM00
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
            
            if ~isscalar(z), error('z must be scalar'); end
            
            % --- Beam Parameters at z ---
            w_z = obj.beam_waist(z);
            R_z = obj.radius_of_curvature(z);
            psi_z = obj.gouy_phase(z);
            
            % --- Amplitude Factor ---
            % E = sqrt(P) * sqrt(2/pi)/w(z) * ...
            % We calculate C so that integral |E|^2 = P_tx
            power_scale = 0;
            if args.P_tx_watts > 0
                power_scale = sqrt(args.P_tx_watts);
            end
            amplitude_factor = (sqrt(2/pi) / w_z) * power_scale;
            
            % --- Radial Profile (Gaussian) ---
            radial_term = exp( -r.^2 ./ w_z^2 );
            
            % --- Phase Terms ---
            % 1. Plane wave propagation
            prop_phase = exp(1i * obj.k * z);
            
            % 2. Gouy Phase shift
            gouy_term = exp(-1i * psi_z);
            
            % 3. Wavefront Curvature
            if isinf(R_z)
                curv_term = 1.0;
            else
                curv_term = exp( -1i * obj.k * r.^2 / (2 * R_z) );
            end
            
            % --- Beam Steering (Tilt) ---
            % Convert r,phi to x,y just for tilt calculation
            x = r .* cos(phi);
            y = r .* sin(phi);
            steering_phase = obj.k * (x * args.beam_tilt_x_rad + y * args.beam_tilt_y_rad);
            
            % --- Noise Injection ---
            phase_noise = 0.0;
            if ~isempty(args.laser_linewidth_kHz) && args.laser_linewidth_kHz > 0
                if isempty(args.symbol_time_s), error('Need symbol_time_s for noise'); end
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
            
            total_noise = exp(1i * (phase_noise + timing_jitter_phase));

            % --- Combine ---
            field = amplitude_factor * ...
                    radial_term .* ...
                    curv_term .* ...
                    gouy_term .* ...
                    prop_phase .* ...
                    exp(1i * steering_phase) .* ...
                    total_noise;
            
            % --- Aperture ---
            if ~isempty(args.tx_aperture_radius)
                mask = (r <= args.tx_aperture_radius);
                field = field .* mask;
            end
        end
        
        function intensity = calculate_intensity(obj, r, phi, z, varargin)
            field = obj.generate_beam_field(r, phi, z, varargin{:});
            intensity = abs(field).^2;
        end
        
        function loss = calculate_coupling_loss(obj, r_aper)
            % Analytic calculation of power loss through a circular aperture
            % Transmission T = 1 - exp(-2 * r_aper^2 / w^2)
            % This assumes the aperture is at the waist for simplicity, 
            % or one must pass the specific w(z).
            % Here we define it generally for the current w0.
            loss_fraction = exp(-2 * r_aper^2 / obj.w0^2);
            loss = -10 * log10(1 - loss_fraction);
        end

        function summary = get_tx_parameters_summary(obj, varargin)
             % Similar utility for consistency
             p_parser = inputParser;
             addParameter(p_parser, 'P_tx_watts', 1.0);
             % (Add other params as needed)
             parse(p_parser, varargin{:});
             summary = p_parser.Results;
             summary.type = 'Gaussian TEM00';
             summary.divergence_mrad = obj.theta_div * 1e3;
        end
    end
end