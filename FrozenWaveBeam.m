classdef FrozenWaveBeam
    % FrozenWaveBeam: Constructs a beam with an arbitrary longitudinal intensity profile.
    % Method: Superposition of co-propagating Bessel beams (Fourier-Bessel expansion).
    % Reference: T. A. Vieira et al., "Frozen waves: experimental realization," 
    % Opt. Lett. 37, 4970 (2012).
    
    properties
        wavelength      % Wavelength
        w0_aperture     % Physical aperture radius (Gaussian apodization)
        L_pattern       % Length of the longitudinal pattern segment
        Q_factor        % Offset parameter (controls transverse spot size)
        N_modes         % Number of Bessel beams to superpose (e.g., 20)
        
        % Internal
        k               % Wavenumber
        A_coeffs        % Calculated weights for each Bessel beam
        kz_n            % Longitudinal wavenumbers
        kr_n            % Transverse wavenumbers
    end
    
    methods
        function obj = FrozenWaveBeam(wavelength, w0_aperture, L_pattern, Q_factor, N_modes, target_func)
            % Constructor
            % L_pattern:   The length over which we define the pattern (e.g., 0.05 m)
            % Q_factor:    Controls the carrier spatial freq (higher = thinner beam)
            % N_modes:     How many harmonics to use (resolution of the pattern)
            % target_func: A function handle @(z) describing desired intensity on-axis
            
            obj.wavelength = wavelength;
            obj.w0_aperture = w0_aperture;
            obj.L_pattern = L_pattern;
            obj.Q_factor = Q_factor;
            obj.N_modes = N_modes;
            obj.k = 2 * pi / wavelength;
            
            % --- 1. Solve for Wavenumbers ---
            % The longitudinal wavevector for the n-th Bessel beam is:
            % kz_n = Q + 2*pi*n / L
            
            n_indices = -N_modes : N_modes;
            obj.kz_n = obj.Q_factor + (2 * pi * n_indices) / obj.L_pattern;
            
            % Verify we aren't evanescent (kz must be < k)
            if any(abs(obj.kz_n) >= obj.k)
                warning('Some modes are evanescent (kz > k). Adjust Q or L.');
                % Filter valid modes
                mask = abs(obj.kz_n) < obj.k;
                obj.kz_n = obj.kz_n(mask);
                n_indices = n_indices(mask);
            end
            
            % Calculate corresponding radial wavenumbers
            % kr = sqrt(k^2 - kz^2)
            obj.kr_n = sqrt(obj.k^2 - obj.kz_n.^2);
            
            % --- 2. Solve for Coefficients A_n ---
            % We want |Sum A_n exp(i kz_n z)|^2 ~ Target(z)
            % This implies A_n are the Fourier Series coefficients of sqrt(Target(z)).
            
            % Numerical integration to find Fourier coefficients
            % F(z) is the complex envelope we want. Let's assume we want real amplitude sqrt(Target).
            
            obj.A_coeffs = zeros(size(n_indices));
            
            % Discretize z for integration
            z_int = linspace(0, L_pattern, 2000);
            dz = z_int(2) - z_int(1);
            
            % Evaluate target shape
            target_amp = sqrt(target_func(z_int));
            
            % Calculate Fourier Coefficients
            % A_n = (1/L) * Integral [ Target(z) * exp(-i * 2*pi*n*z/L) ] dz
            for i = 1:length(n_indices)
                n = n_indices(i);
                basis = exp(-1i * 2 * pi * n * z_int / obj.L_pattern);
                integrand = target_amp .* basis;
                obj.A_coeffs(i) = (1 / obj.L_pattern) * sum(integrand) * dz;
            end
        end
        
        function field = generate_beam_field(obj, r, phi, z, varargin)
            p_parser = inputParser;
            addRequired(p_parser, 'r');
            addRequired(p_parser, 'phi');
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            parse(p_parser, r, phi, z, varargin{:});
            args = p_parser.Results;
            
            if ~isscalar(z), error('z must be scalar'); end
            
            % Superposition of Bessel Beams
            % E(r,z) = Sum [ A_n * J0(kr_n * r) * exp(i * kz_n * z) ]
            
            total_field = zeros(size(r));
            
            % Since this simulates a physical experiment, we must apply 
            % the Gaussian aperture to ALL Bessel components.
            % Ideally, J0 should be J0_approx (Bessel-Gauss).
            % We apply the envelope at the end for simplicity, valid for paraxial.
            
            for i = 1:length(obj.kz_n)
                An = obj.A_coeffs(i);
                kr = obj.kr_n(i);
                kz = obj.kz_n(i);
                
                % Bessel term
                bessel_term = besselj(0, kr * r);
                
                % Prop term
                prop_term = exp(1i * kz * z);
                
                total_field = total_field + An * bessel_term * prop_term;
            end
            
            % Apply physical aperture (Gaussian Apodization)
            % This limits the distance the frozen wave persists, but makes it physical.
            gauss_env = exp( -r.^2 / obj.w0_aperture^2 );
            
            scale = 1.0;
            if args.P_tx_watts > 0, scale = sqrt(args.P_tx_watts); end
            
            field = scale * total_field .* gauss_env;
        end
        
        function intensity = calculate_intensity(obj, r, phi, z, varargin)
            field = obj.generate_beam_field(r, phi, z, varargin{:});
            intensity = abs(field).^2;
        end
    end
end