classdef FlatTopBeam
    % FlatTopBeam: Models a beam with uniform intensity (Top-Hat).
    % Implementation: Gori's Flattened Gaussian Beam (FGB) expansion.
    % Reference: F. Gori, "Flattened gaussian beams," Opt. Commun. 107, 335 (1994).
    
    properties
        N_order         % Flatness order (N=0 is Gaussian, N=20 is sharp Top-Hat)
        w0              % Waist parameter (Controls the width of the Top-Hat)
        wavelength      % Wavelength
        
        % Derived
        k               % Wavenumber
        z0              % Fundamental Rayleigh range
        cn_coeffs       % Expansion coefficients for the LG superposition
    end
    
    methods
        function obj = FlatTopBeam(N_order, w0, wavelength)
            % Constructor
            if N_order < 0, error('Order N must be >= 0'); end
            
            obj.N_order = N_order;
            obj.w0 = w0;
            obj.wavelength = wavelength;
            
            obj.k = 2 * pi / wavelength;
            obj.z0 = (pi * w0^2) / wavelength;
            
            % --- Precompute Gori Coefficients ---
            % c_n = (-1)^n * Sum_{m=n}^N [ Binomial(m,n) * (1/2^m) ]
            % These weights determine how much of each LG mode is needed.
            
            obj.cn_coeffs = zeros(1, N_order + 1);
            for n = 0:N_order
                sum_val = 0;
                for m = n:N_order
                    sum_val = sum_val + nchoosek(m, n) * (1/2^m);
                end
                obj.cn_coeffs(n+1) = (-1)^n * sum_val;
            end
        end
        
        function field = generate_beam_field(obj, r, phi, z, varargin)
            p_parser = inputParser;
            addRequired(p_parser, 'r');
            addRequired(p_parser, 'phi');
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            addParameter(p_parser, 'tx_aperture_radius', []);
            parse(p_parser, r, phi, z, varargin{:});
            args = p_parser.Results;
            
            if ~isscalar(z), error('z must be scalar'); end
            
            % --- Superposition Logic ---
            % E(r,z) = Sum_{n=0}^N [ c_n * LG_n^0(r,z) ]
            % We generate the base LG_n^0 modes analytically here for speed.
            
            total_field = zeros(size(r));
            
            % Common parameters for all modes
            % w(z), R(z) are same for all LG modes sharing same w0
            w_z = obj.w0 * sqrt(1 + (z/obj.z0)^2);
            R_z = z * (1 + (obj.z0/z)^2);
            if abs(z) < 1e-12, R_z = inf; end
            
            % Base Gaussian Envelope
            % Exp( -r^2/w^2 ) * Exp( -ik r^2 / 2R ) * Exp( ikz )
            r_sq = r.^2;
            gauss_part = (obj.w0 / w_z) .* exp( -r_sq ./ w_z^2 );
            
            if isinf(R_z)
                phase_part = exp(1i * obj.k * z);
            else
                phase_part = exp(1i * obj.k * z) .* exp( -1i * obj.k * r_sq / (2 * R_z) );
            end
            
            base_field = gauss_part .* phase_part;
            
            % Laguerre Polynomial Argument
            % arg = 2 * r^2 / w(z)^2
            arg_L = 2 * r_sq ./ w_z^2;
            
            for n = 0:obj.N_order
                coeff = obj.cn_coeffs(n+1);
                
                % 1. Gouy Phase for mode n (l=0)
                % Phase = (2n + 1) * atan(z/z0)
                gouy_phase = exp( -1i * (2*n + 1) * atan(z / obj.z0) );
                
                % 2. Laguerre Polynomial L_n^0(x)
                % Using static numerical solver for efficiency
                L_val = FlatTopBeam.numerical_laguerre(n, arg_L);
                
                % Summation
                total_field = total_field + coeff .* L_val .* gouy_phase;
            end
            
            % Apply Base Field
            total_field = total_field .* base_field;
            
            % Power Scaling
            % Note: FGB normalization is tricky. 
            % At N=infinity, Intensity = 1. At N=0, Intensity = 1.
            % We apply simple P_tx scaling assuming the user wants peak-like scaling
            % or use a numerical norm factor if strict power conservation is needed.
            % Here we maintain the relative shape defined by Gori.
            
            scale = 1.0;
            if args.P_tx_watts > 0, scale = sqrt(args.P_tx_watts); end
            
            field = scale * total_field;
            
            if ~isempty(args.tx_aperture_radius)
                field = field .* (r <= args.tx_aperture_radius);
            end
        end
        
        function intensity = calculate_intensity(obj, r, phi, z, varargin)
            field = obj.generate_beam_field(r, phi, z, varargin{:});
            intensity = abs(field).^2;
        end
    end
    
    methods (Static)
        function L_val = numerical_laguerre(n, x)
            % Fast L_n^0(x) computation via recurrence
            if n == 0, L_val = ones(size(x)); return; end
            if n == 1, L_val = 1 - x; return; end
            
            L_km1 = ones(size(x));
            L_k = 1 - x;
            
            for k = 1:(n-1)
                % (k+1) L_{k+1} = (2k + 1 - x)L_k - k L_{k-1}
                L_next = ( (2*k + 1 - x).*L_k - k*L_km1 ) / (k+1);
                L_km1 = L_k;
                L_k = L_next;
            end
            L_val = L_k;
        end
    end
end