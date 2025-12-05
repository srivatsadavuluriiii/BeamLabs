% run_gaussian_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 1550e-9;   % Telecom C-band
W0 = 2.0e-3;            % 2 mm beam waist (collimated output)
GRID_SIZE = 300;
MAX_RADIUS_MM = 6.0;    % View window
SAVE_FIGURE = false;
PLOT_DIR = 'plots_gaussian_matlab';

if ~exist(PLOT_DIR, 'dir') && SAVE_FIGURE; mkdir(PLOT_DIR); end

fprintf('\n--- Fundamental Gaussian Beam (TEM00) Analysis ---\n');

% Instantiate
beam = GaussianBeam(WAVELENGTH, W0);

fprintf('Beam Parameters:\n');
fprintf('  Wavelength: %.0f nm\n', beam.wavelength*1e9);
fprintf('  Waist w0: %.3f mm\n', beam.w0*1e3);
fprintf('  Rayleigh Range z0: %.2f m\n', beam.z0);
fprintf('  Far-Field Divergence: %.2f mrad\n', beam.theta_div*1e3);

%% Visual Setup
r_max_m = MAX_RADIUS_MM * 1e-3;
x_vec = linspace(-r_max_m, r_max_m, GRID_SIZE);
y_vec = linspace(-r_max_m, r_max_m, GRID_SIZE);
[X, Y] = meshgrid(x_vec, y_vec);
[PHI, R_grid] = cart2pol(X, Y);

%% 1. Near Field (Waist) vs Far Field
z_near = 0;         % At waist
z_far = 10 * beam.z0; % Deep far field

% Calculate Fields
field_near = beam.generate_beam_field(R_grid, PHI, z_near);
field_far = beam.generate_beam_field(R_grid, PHI, z_far);

% Normalize intensities for visualization (since far field is much dimmer)
int_near = abs(field_near).^2;
int_far = abs(field_far).^2; 
int_far_norm = int_far / max(int_far(:)); % Normalize peak to 1

% Calculate Phases
phase_near = angle(field_near);
phase_far = angle(field_far);

%% Plotting

% --- Figure 1: Intensity & Phase Evolution ---
figure('Position', [100, 100, 1000, 700], 'Color', 'w');

% Near Field Intensity
subplot(2, 2, 1);
imagesc(x_vec*1e3, y_vec*1e3, int_near);
axis xy equal tight; colormap('hot'); colorbar;
title('Intensity @ z=0 (Waist)');
xlabel('x [mm]'); ylabel('y [mm]');

% Near Field Phase (Should be Flat)
subplot(2, 2, 2);
imagesc(x_vec*1e3, y_vec*1e3, phase_near);
axis xy equal tight; 
try colormap(gca, twilight); catch, colormap(gca, jet); end 
colorbar;
title('Phase @ z=0 (Planar Wavefront)');

% Far Field Intensity
subplot(2, 2, 3);
imagesc(x_vec*1e3, y_vec*1e3, int_far_norm);
axis xy equal tight; colormap('hot'); colorbar;
title(sprintf('Norm. Intensity @ z=%.1f m (Expanded)', z_far));

% Far Field Phase (Should be Wrapped Spherical Rings)
subplot(2, 2, 4);
imagesc(x_vec*1e3, y_vec*1e3, phase_far);
axis xy equal tight; 
try colormap(gca, twilight); catch, colormap(gca, jet); end 
colorbar;
title('Phase @ Far Field (Spherical Wavefront)');


% --- Figure 2: Beam Caustic (Side View) ---
% This verifies the divergence angle
z_prop = linspace(0, 2*beam.z0, 200);
r_prop = linspace(-3*beam.w0, 3*beam.w0, 200);
[Z_mesh, R_mesh] = meshgrid(z_prop, r_prop);

% Vectorized calculation for plot speed
% (Manual calc since generate_beam_field enforces scalar z for rigorous phase noise)
w_z_vec = beam.w0 * sqrt(1 + (Z_mesh./beam.z0).^2);
int_prop = (beam.w0 ./ w_z_vec).^2 .* exp( -2 * R_mesh.^2 ./ w_z_vec.^2 );

figure('Position', [150, 150, 800, 500], 'Color', 'w');
imagesc(z_prop, r_prop*1e3, int_prop);
colormap('hot');
axis xy;
xlabel('Propagation Distance z [m]');
ylabel('Radial Position r [mm]');
title('Gaussian Beam Caustic (Side View)');
hold on;

% Plot the 1/e^2 radius lines
w_analytical = beam.beam_waist(z_prop);
plot(z_prop, w_analytical*1e3, 'c--', 'LineWidth', 2);
plot(z_prop, -w_analytical*1e3, 'c--', 'LineWidth', 2);
legend('Intensity', 'Beam Radius w(z)');

% Mark Rayleigh Range
xline(beam.z0, 'g:', 'Rayleigh Range z_R');

fprintf('Analysis Complete.\n');