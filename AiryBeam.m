classdef AiryBeam
    % AiryBeam: Models the 1D or 2D Finite-Energy Airy Beam.
    % Based on Siviloglou & Christodoulides, Opt. Lett. 32, 979 (2007).
    
    properties
        x0              % Transverse scale (main lobe width parameter)
        a               % Apodization (decay) parameter (usually 0 < a << 1)
        wavelength      % Wavelength in meters
        
        % Derived Constants
        k               % Wavenumber
        z0_airy         % Characteristic diffraction length (different from Rayleigh range!)
    end
    
    methods
        function obj = AiryBeam(wavelength, x0, a)
            % Constructor
            % x0: Main lobe scale (e.g., 10-100 microns). Controls bending rate.
            % a:  Truncation factor (e.g., 0.05). a=0 is ideal (infinite energy).
            
            if a < 0, error('Apodization parameter "a" must be >= 0'); end
            
            obj.wavelength = wavelength;
            obj.x0 = x0;
            obj.a = a;
            
            obj.k = 2 * pi / wavelength;
            
            % Characteristic diffraction length for Airy beams:
            % z0 = 2 * k * x0^2
            obj.z0_airy = 2 * obj.k * obj.x0^2; 
        end
        
        function x_deflection = trajectory(obj, z)
            % Calculates the transverse deflection of the main lobe
            % x_d = z^2 / (4 * k^2 * x0^3)
            % Note: This is the parabolic path.
            
            % Normalized distance xi = z / (k x0^2)
            xi = z / (obj.k * obj.x0^2);
            
            % Deflection in normalized units is xi^2 / 4
            % Convert back to meters: x = x0 * (xi^2 / 4)
            x_deflection = obj.x0 * (xi.^2) / 4;
        end
        
        function field = generate_beam_field(obj, x, y, z, varargin)
            % Generates 2D Airy beam E(x,y,z) = Ai(x,z) * Ai(y,z)
            % Parabolic acceleration occurs in both x and y directions (usually).
            % Use 'is_1d' flag if you only want a 1D sheet (invariant in y).
            
            p_parser = inputParser;
            addRequired(p_parser, 'x');
            addRequired(p_parser, 'y');
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            addParameter(p_parser, 'is_1d', false); % If true, ignores y-variation (light sheet)
            addParameter(p_parser, 'tx_aperture_radius', []);
            
            parse(p_parser, x, y, z, varargin{:});
            args = p_parser.Results;
            
            if ~isscalar(z), error('z must be scalar'); end
            
            % --- Normalized Units ---
            % sX = x / x0
            % sY = y / x0
            % xi = z / (k * x0^2)
            
            sX = x / obj.x0;
            sY = y / obj.x0;
            xi = z / (obj.k * obj.x0^2);
            
            % --- 1D Airy Function Calculation (function internal below) ---
            % Field_X = Ai( sX - (xi/2)^2 + i*a*xi ) * ...
            %           exp( i*sX*xi/2 - i*xi^3/12 + ... )
            
            u_x = obj.calculate_1d_airy(sX, xi);
            
            if args.is_1d
                u_y = ones(size(y)); % Constant in Y
                % Normalize for 1D sheet power? keeping simple for now.
            else
                u_y = obj.calculate_1d_airy(sY, xi);
            end
            
            % Total Field
            field_norm = u_x .* u_y;
            
            % Power Scaling (Approximate)
            % Airy beam power normalization is complex because it depends on 'a'.
            % We apply a simple sqrt(P) scalar.
            scale = 1.0;
            if args.P_tx_watts > 0
                scale = sqrt(args.P_tx_watts); 
            end
            
            field = scale * field_norm;
            
            % Aperture
            if ~isempty(args.tx_aperture_radius)
                r = sqrt(x.^2 + y.^2);
                mask = (r <= args.tx_aperture_radius);
                field = field .* mask;
            end
        end
        
        function intensity = calculate_intensity(obj, x, y, z, varargin)
            field = obj.generate_beam_field(x, y, z, varargin{:});
            intensity = abs(field).^2;
        end
        
        function [u_1d] = calculate_1d_airy(obj, s, xi)
            % Implements the analytical Finite-Energy Airy solution
            % Ref: Siviloglou Eq (4)
            % s: Normalized transverse coordinate (x/x0)
            % xi: Normalized propagation distance
            % a: Decay parameter
            
            a_param = obj.a;
            
            % Argument for the Airy function:
            % arg = s - (xi/2)^2 + i * a * xi
            airy_arg = s - (xi/2)^2 + 1i * a_param * xi;
            
            % Argument for the Exponential:
            % term1 = i * s * xi / 2
            % term2 = - i * xi^3 / 12
            % term3 = i * a^2 * xi / 2
            % term4 = - a * s
            % term5 = a * xi^2 / 2
            
            exp_arg = (1i * s .* xi / 2) ...
                    - (1i * xi^3 / 12) ...
                    + (1i * a_param^2 * xi / 2) ...
                    - (a_param * s) ...
                    + (a_param * xi^2 / 2);
            
            % Calculate
            % Note: MATLAB's airy(0, z) computes Ai(z)
            val_Ai = airy(0, airy_arg);
            
            u_1d = val_Ai .* exp(exp_arg);
        end
        
        function summary = get_beam_summary(obj)
            summary.type = 'Finite-Energy Airy Beam';
            summary.x0 = obj.x0;
            summary.a = obj.a;
            summary.z0_airy = obj.z0_airy;
            % Calculate deflection at 1 diffraction length
            summary.deflection_at_z0 = obj.trajectory(obj.z0_airy);
        end
    end
end