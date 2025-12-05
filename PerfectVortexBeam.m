classdef PerfectVortexBeam
    % PerfectVortexBeam: Models a beam with OAM but fixed ring radius.
    % Implementation based on the Fourier Transform of a Bessel-Gaussian beam.
    % Reference: Vaity & Rusch, Opt. Lett. 40, 597 (2015).
    
    properties
        l               % Topological charge (OAM order)
        r_ring          % Desired fixed radius of the ring (meters)
        w0              % Width of the ring (approximate Gaussian thickness)
        wavelength      % Wavelength in meters
        
        % Derived
        k               % Wavenumber
        f_lens          % Virtual focal length (used for numerical generation context)
    end
    
    methods
        function obj = PerfectVortexBeam(l, r_ring, w0, wavelength)
            % Constructor
            % l:      OAM order (integer)
            % r_ring: The fixed radius of the bright ring
            % w0:     The thickness of the ring (waist of the generating Gaussian)
            
            obj.l = l;
            obj.r_ring = r_ring;
            obj.w0 = w0;
            obj.wavelength = wavelength;
            obj.k = 2 * pi / wavelength;
            
            % We define an arbitrary focal length for the "Fourier Transform" 
            % propagation model, usually normalized to 1m or similar for calc.
            obj.f_lens = 1.0; 
        end
        
        function field = generate_beam_field(obj, r, phi, z, varargin)
            % Generates the POV field.
            % Note: Analytical POV expressions involve complex Modified Bessel 
            % functions I_l. We implement the rigorous analytical form.
            % E(r) ~ i^(l-1) * (w0/w(z)) * exp(...) * I_l(...)
            
            p_parser = inputParser;
            addRequired(p_parser, 'r');
            addRequired(p_parser, 'phi');
            addRequired(p_parser, 'z');
            addParameter(p_parser, 'P_tx_watts', 1.0);
            addParameter(p_parser, 'tx_aperture_radius', []);
            
            parse(p_parser, r, phi, z, varargin{:});
            args = p_parser.Results;
            
            if ~isscalar(z), error('z must be scalar'); end
            
            % --- Physics Constants ---
            % The "Perfect" Vortex is formed at the focus of a lens acting on a Bessel beam.
            % Here we model the field *at* that focal plane (z=0 for our class)
            % and its propagation away from it.
            
            % Standard Gaussian waist expansion logic applies to the "thickness" of the ring
            z0_ring = (pi * obj.w0^2) / obj.wavelength;
            w_z = obj.w0 * sqrt(1 + (z / z0_ring)^2);
            
            % Complex q-parameter stuff for the envelope
            % This is an approximation of the POV analytical solution:
            % E(r,z) = A * exp(-r^2/wz^2) * exp(-Rr^2/wz^2) * I_n(2*r*Rr/wz^2)
            
            % 1. Amplitude Scaling
            scale_factor = 1 / w_z; % Amplitude drops as ring thickens
            
            % 2. Gaussian Envelopes
            % The beam is fundamentally a Gaussian ring at r_ring
            arg_gauss = - (r.^2 + obj.r_ring^2) ./ (w_z^2);
            val_gauss = exp(arg_gauss);
            
            % 3. Modified Bessel Function of the First Kind (I_l)
            % This is the core engine of the POV.
            % Note: For large arguments, I_l(x) grows exponentially.
            % However, val_gauss decays exponentially.
            % We must use scaled Bessel to avoid Infinity * 0 NaN errors.
            % besseli(nu, Z, 1) computes besseli(nu, Z) .* exp(-abs(real(Z)))
            
            arg_bessel = (2 * r * obj.r_ring) ./ (w_z^2);
            
            % Use scaled Bessel I to maintain numerical stability
            % I_l(x) * exp(-x) * exp(x) -> we use the scaled version
            % Real argument here, so exp(-abs(real(x))) = exp(-x)
            val_bessel_scaled = besseli(obj.l, arg_bessel, 1); 
            
            % Re-apply the exponential we "borrowed" for stability
            % The total exponential term becomes: 
            % exp( -(r^2 + R^2)/w^2 ) * exp( 2*r*R/w^2 ) 
            % = exp( -(r - R)^2 / w^2 )  <-- A Gaussian centered at R!
            
            % So mathematically, we construct the field using the stable form:
            stable_envelope = exp( - (r - obj.r_ring).^2 ./ w_z^2 );
            
            % 4. Phase Terms
            azimuthal_phase = exp(1i * obj.l * phi);
            gouy_phase = exp(-1i * atan(z / z0_ring)); % Ring Gouy phase
            prop_phase = exp(1i * obj.k * z);
            
            % Curvature (applies to the ring width)
            if abs(z) > 1e-12
                R_curv = z * (1 + (z0_ring/z)^2);
                curv_phase = exp(-1i * obj.k * (r - obj.r_ring).^2 / (2 * R_curv));
            else
                curv_phase = 1.0;
            end

            % Power Normalization
            % Rough approx: 2*pi*R * w0 is effective area
            % We apply user power
            p_scale = 1.0;
            if args.P_tx_watts > 0, p_scale = sqrt(args.P_tx_watts); end
            
            % --- Assembly ---
            % Field ~ Gaussian(r-R) * ScaledBessel * Phase
            % Note: For r ~ R, the scaled Bessel function approaches a constant 1/sqrt(2*pi*arg).
            % We ignore the amplitude variation of the Bessel for pure visualization
            % as the Gaussian envelope dominates.
            
            field = p_scale * scale_factor * ...
                    stable_envelope .* ...
                    val_bessel_scaled .* ...
                    azimuthal_phase .* ...
                    gouy_phase .* ...
                    curv_phase .* ...
                    prop_phase;
            
            if ~isempty(args.tx_aperture_radius)
                field = field .* (r <= args.tx_aperture_radius);
            end
        end
        
        function intensity = calculate_intensity(obj, r, phi, z, varargin)
            field = obj.generate_beam_field(r, phi, z, varargin{:});
            intensity = abs(field).^2;
        end
    end
end