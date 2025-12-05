classdef HelicoConicalBeam
    % HelicoConicalBeam: Models a beam with spiral phase-radial coupling.
    % Reference: C. A. Alonzo et al., "Helico-conical optical beams," 
    % Opt. Lett. 30, 3050 (2005).
    
    properties
        l               % Topological charge (OAM) - Controls spiral head
        K               % Conical parameter - Controls spiral tail winding
        r0              % Radial normalization scale
        wavelength      % Wavelength
        w0_gauss        % Gaussian apodization radius
        
        % Derived
        k               % Wavenumber
    end
    
    methods
        function obj = HelicoConicalBeam(l, K, r0, wavelength, w0_gauss)
            % Constructor
            obj.l = l;
            obj.K = K;
            obj.r0 = r0;
            obj.wavelength = wavelength;
            obj.w0_gauss = w0_gauss;
            obj.k = 2 * pi / wavelength;
        end
        
        function field = generate_source_field(obj, r, phi, varargin)
            % Generates the field at the source plane (z=0).
            % Phase = l*phi + K*(r/r0 - 1)
            
            p_parser = inputParser;
            addRequired(p_parser, 'r');
            addRequired(p_parser, 'phi');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            parse(p_parser, r, phi, varargin{:});
            args = p_parser.Results;
            
            % --- 1. Helico-Conical Phase ---
            % We normalize r/r0.
            % Note: The term '-1' is a constant phase offset often included 
            % in literature to set the zero-phase point at r=r0.
            phase_helico = obj.l * phi + obj.K * (r ./ obj.r0 - 1);
            
            % --- 2. Gaussian Apodization ---
            % Crucial for physical realization and FFT stability (prevents hard edge ringing)
            amp_gauss = exp( - (r.^2) ./ obj.w0_gauss^2 );
            
            % --- 3. Assembly ---
            scale_factor = 1.0;
            if args.P_tx_watts > 0
                scale_factor = sqrt(args.P_tx_watts);
            end
            
            field = scale_factor * amp_gauss .* exp(1i * phase_helico);
        end
        
        function intensity = calculate_intensity(obj, r, phi, varargin)
            field = obj.generate_source_field(r, phi, varargin{:});
            intensity = abs(field).^2;
        end
    end
    
    methods (Static)
        function field_out = propagate_to_focus(field_in)
            % Performs the Optical Fourier Transform (Lens focusing).
            % E_focus(u, v) ~ FFT2(E_source(x, y))
            %
            % RIGOROUS SHIFTING SEQUENCE:
            % 1. ifftshift: Moves zero-frequency (center of grid) to index (1,1)
            % 2. fft2:      Computes FFT
            % 3. fftshift:  Moves DC component back to the center of the array
            
            field_out = fftshift(fft2(ifftshift(field_in)));
            
            % Note on Amplitude Scaling:
            % FFT conserves Parseval's energy sum, but pixel size changes.
            % For pure shape visualization, raw FFT output is sufficient.
            % For physical units, one would multiply by (dx*dy / (i*lambda*f)).
        end
    end
end