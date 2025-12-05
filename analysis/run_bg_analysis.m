% run_bg_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 1550e-9;  
W0 = 2.5e-3;         % Gaussian Envelope Waist (2.5 mm)
                     % Note: Usually w0 is large to support a long Bessel zone
L_MODE = 2;          % Topological charge (0 = central spot, >0 = dark core)

% Defining BETA (The transverse wavenumber)
% High beta = tight rings, short diffraction-free range
% Low beta = wide rings, long diffraction-free range
% Let's define it to have a Core Radius of ~100 microns
target_core_radius = 100e-6; 
BETA = 2.4048 / target_core_radius; 

GRID_SIZE = 250;
MAX_RADIUS_MM = 4.0; % View window
SAVE_FIGURE = false; % Set true to save
PLOT_DIR = 'plots_bg_matlab';

if ~exist(PLOT_DIR, 'dir') && SAVE_FIGURE
    mkdir(PLOT_DIR);
end

fprintf('\n--- Bessel-Gaussian (BG) Beam Analysis ---\n');

% Instantiate Beam
beam = BesselGaussianBeam(L_MODE, BETA, WAVELENGTH, W0);

fprintf('Beam Parameters:\n');
fprintf('  Lambda: %.0f nm\n', beam.wavelength*1e9);
fprintf('  Gaussian Waist (w0): %.2f mm\n', beam.w0*1e3);
fprintf('  Bessel Order (l): %d\n', beam.l);
fprintf('  Transverse K (beta): %.1f rad/m\n', beam.beta);
fprintf('  Geometric Cone Angle: %.3f deg\n', rad2deg(beam.alpha_cone));
fprintf('  Diffraction-Free Range (z_max): %.2f m\n', beam.z_max);

%% Propagation Check
% We want to check points inside and outside the "Bessel Zone"
z_list = [0, beam.z_max*0.5, beam.z_max, beam.z_max*1.5, beam.z_max*3.0];
beam.propagation_summary(z_list);

%% Visual Analysis Setup
r_max_m = MAX_RADIUS_MM * 1e-3;
x_vec = linspace(-r_max_m, r_max_m, GRID_SIZE);
y_vec = linspace(-r_max_m, r_max_m, GRID_SIZE);
[X, Y] = meshgrid(x_vec, y_vec);
[PHI, R_grid] = cart2pol(X, Y);

% Field at z = 0.5 * z_max (Deep in diffraction-free zone)
zField = beam.z_max * 0.5;
field_z = beam.generate_beam_field(R_grid, PHI, zField);
intensity_z = abs(field_z).^2;
phase_z = angle(field_z);

% Radial Profile
r_range = linspace(0, r_max_m, GRID_SIZE*2);
intensity_profile = beam.calculate_intensity(r_range, zeros(size(r_range)), zField);

%% Plotting

% --- Figure 1: Transverse Intensity (Central Spot/Ring) ---
figure('Position', [100, 100, 900, 600], 'Color', 'w');
subplot(1,2,1);
imagesc(x_vec*1e3, y_vec*1e3, intensity_z);
axis xy equal tight;
colormap('hot');
title(sprintf('Intensity @ z=%.1f m (In Bessel Zone)', zField));
xlabel('x [mm]'); ylabel('y [mm]');
colorbar;

subplot(1,2,2);
plot(r_range*1e3, intensity_profile, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Radius [mm]'); ylabel('Intensity');
title('Radial Profile');
xlim([0, MAX_RADIUS_MM]);

% --- Figure 2: Phase (Check OAM) ---
figure('Position', [150, 150, 600, 500], 'Color', 'w');
imagesc(x_vec*1e3, y_vec*1e3, phase_z);
axis xy equal tight;
% Custom twilight map logic
try colormap(twilight); catch, colormap(hsv); end 
colorbar;
title(sprintf('Phase @ z=%.1f m (l=%d)', zField, beam.l));
xlabel('x [mm]'); ylabel('y [mm]');

% --- Figure 3: Longitudinal Propagation (The "Triangle" of Existence) ---
% This is the most important plot for Bessel beams.
% We plot from z=0 past z_max to see the core dissolve.

z_max_plot = beam.z_max * 2.0;
num_z_steps = 200;
z_array = linspace(0, z_max_plot, num_z_steps);

% We focus strictly on the center to see the core evolution
r_max_long = 0.5e-3; % 0.5 mm radius view
num_r_steps = 200;
r_array_long = linspace(-r_max_long, r_max_long, num_r_steps);

intensity_long = zeros(num_r_steps, num_z_steps);

fprintf('\nCalculating Longitudinal Map...\n');
for i = 1:num_z_steps
    z_val = z_array(i);
    % Calculate 1D slice at x=0 (phi=0 or pi)
    % Pass phi=0 for positive r, phi=pi for negative r if strictly polar, 
    % but since r can be treated as coordinate in calculate_intensity:
    intensity_long(:, i) = beam.calculate_intensity(abs(r_array_long), zeros(size(r_array_long)), z_val);
end

figure('Position', [200, 200, 1000, 600], 'Color', 'w');
imagesc(z_array, r_array_long*1e3, intensity_long);
axis xy;
colormap('hot');
colorbar;
xlabel('Propagation Distance z [m]');
ylabel('Radial Position r [mm]');
title(sprintf('Bessel-Gauss Propagation (z_{max} = %.2f m)', beam.z_max));

% Draw the theoretical geometrical shadow lines
% The core exists where |r| < w0 - z * tan(alpha)
% Or more simply, the cone edge: r_edge = w0 - z * (beta/k)
slope = beam.beta / beam.k;
upper_line = (beam.w0 - z_array * slope) * 1e3;
lower_line = -(beam.w0 - z_array * slope) * 1e3;

hold on;
xline(beam.z_max, 'g--', 'Max Range (Geom)', 'LineWidth', 2);
yline(0, 'w:', 'LineWidth', 0.5);

% Visualize the "Cone of Silence" (Shadow zone)
% Plotting these lines helps visualize why the beam dies
plot(z_array, upper_line, 'c--', 'LineWidth', 1);
plot(z_array, lower_line, 'c--', 'LineWidth', 1);

ylim([-r_max_long*1e3, r_max_long*1e3]);

fprintf('Analysis Complete.\n');