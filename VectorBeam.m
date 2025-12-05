classdef VectorBeam
    % VectorBeam: Models beams with spatially varying polarization.
    % Supports Radial, Azimuthal, and Generalized Cylindrical Vector Beams (CVBs).
    
    properties
        type            % 'radial', 'azimuthal', or 'hybrid'
        order           % Polarization order P (usually 1 for standard Rad/Azim)
        wavelength      % Wavelength
        w0              % Beam waist
        
        % Derived
        k               % Wavenumber
        z0              % Rayleigh range
    end
    
    methods
        function obj = VectorBeam(type, wavelength, w0, order)
            if nargin < 4, order = 1; end
            
            valid_types = {'radial', 'azimuthal', 'hybrid'};
            if ~ismember(lower(type), valid_types)
                error('Type must be radial, azimuthal, or hybrid.');
            end
            
            obj.type = lower(type);
            obj.wavelength = wavelength;
            obj.w0 = w0;
            obj.order = order;
            
            obj.k = 2 * pi / wavelength;
            obj.z0 = (pi * w0^2) / wavelength;
        end
        
        function [Ex, Ey, Ez] = generate_vector_fields(obj, r, phi, z, varargin)
            % Returns the vector components.
            % Note: We return Ex, Ey (transverse) and Ez (longitudinal).
            % Ez is usually 0 in paraxial approximation, but for Radial beams
            % it is non-zero near focus. We will model the Transverse parts
            % rigorously here.
            
            p_parser = inputParser;
            addRequired(p_parser, 'r');
            addRequired(p_parser, 'phi');
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            parse(p_parser, r, phi, z, varargin{:});
            args = p_parser.Results;
            
            % --- 1. Base Scalar Profile ---
            % CVBs are usually built on a doughnut-shaped LG_01 mode amplitude
            % because the center must be dark (singularity in polarization).
            
            w_z = obj.w0 * sqrt(1 + (z/obj.z0)^2);
            R_z = z * (1 + (obj.z0/z)^2);
            if abs(z) < 1e-12, R_z = inf; end
            psi_z = 2 * atan(z / obj.z0); % Gouy phase for p=0, l=1
            
            % LG_01 Amplitude Envelope (Doughnut)
            % E ~ r * exp(-r^2/w^2)
            amp_factor = sqrt(2/pi) / w_z * (sqrt(2)*r/w_z) .* exp(-r.^2/w_z^2);
            
            % Phase terms
            phase_factor = exp(1i*obj.k*z) .* exp(-1i*psi_z);
            if ~isinf(R_z)
                phase_factor = phase_factor .* exp(-1i*obj.k*r.^2/(2*R_z));
            end
            
            scalar_profile = amp_factor .* phase_factor;
            
            % Power Scaling
            if args.P_tx_watts > 0
                scalar_profile = scalar_profile * sqrt(args.P_tx_watts);
            end
            
            % --- 2. Polarization Matrix ---
            % We construct the Jones Vectors spatially.
            
            P = obj.order;
            
            if strcmp(obj.type, 'radial')
                % Radial: Polarization points along r-hat
                % x-component = cos(phi), y-component = sin(phi)
                pol_x = cos(P * phi);
                pol_y = sin(P * phi);
                
            elseif strcmp(obj.type, 'azimuthal')
                % Azimuthal: Polarization points along phi-hat
                % x-component = -sin(phi), y-component = cos(phi)
                pol_x = -sin(P * phi);
                pol_y = cos(P * phi);
                
            else % Hybrid (Generalized)
                % Arbitrary rotation angle phi_0 (let's say 45 deg)
                phi_0 = pi/4;
                pol_x = cos(P * phi + phi_0);
                pol_y = sin(P * phi + phi_0);
            end
            
            % --- 3. Vector Components ---
            Ex = scalar_profile .* pol_x;
            Ey = scalar_profile .* pol_y;
            
            % Ez is significant only if strongly focused (NA > 0.7).
            % In paraxial regime, Ez ~ 0.
            Ez = zeros(size(Ex)); 
        end
        
        function intensity = calculate_intensity(obj, r, phi, z, varargin)
            [Ex, Ey, ~] = obj.generate_vector_fields(r, phi, z, varargin{:});
            % Total Intensity is sum of component intensities
            intensity = abs(Ex).^2 + abs(Ey).^2;
        end
    end
end