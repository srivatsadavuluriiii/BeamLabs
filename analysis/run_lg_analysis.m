% run_lg_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 1550e-9;
W0 = 25e-3;
P_MODE = 0;
L_MODE = 2;
GRID_SIZE = 200;
MAX_RADIUS_MM = 75;
SAVE_FIGURE = true;
PLOT_DIR = 'plots_matlab';

if ~exist(PLOT_DIR, 'dir') && SAVE_FIGURE
    mkdir(PLOT_DIR);
end


fprintf('\n--- Laguerre Gaussian Beam Analysis ---\n');

% Instantiate Beam
beam = LaguerreGaussianBeam(P_MODE, L_MODE, WAVELENGTH, W0);
fprintf('Beam: LG p=%d l=%d, lambda=%.0f nm, w0=%.1f mm\n', ...
    beam.p, beam.l, beam.wavelength*1e9, beam.w0*1e3);

%% Beam Parameters
[theta_0, theta_eff] = beam.effective_divergence_angle();
fprintf('Ideal Divergence (theta_0): %.2f urad\n', theta_0 * 1e6);
fprintf('Effective Divergence (Theta): %.2f urad\n', theta_eff * 1e6);
fprintf('Fundamental Rayleigh Range (z0): %.2f m\n', beam.z0_fundamental);

params_at_zR = beam.get_beam_parameters(beam.z_R);
fprintf('\nParams at z=z_R (%.1f m):\n', beam.z_R);
fprintf('  w(z_R) = %.2f mm\n', params_at_zR.w_z * 1e3);
fprintf('  R(z_R) = %.2f m\n', params_at_zR.R_z);

z_list = [0, beam.z_R, 2*beam.z_R, 3*beam.z_R, 1000, 1200];
beam.propagation_summary(z_list);

%% Transmitter Config & Summary
P_tx = 2.0;
laser_linewidth = 10.0; % kHz
timing_jitter = 5.0;    % ps
tx_aperture = 50e-3;    % m
tilt_x = deg2rad(0.1);
tilt_y = 0;

tx_summary = beam.get_tx_parameters_summary(...
    'P_tx_watts', P_tx, ...
    'laser_linewidth_kHz', laser_linewidth, ...
    'timing_jitter_ps', timing_jitter, ...
    'tx_aperture_radius', tx_aperture, ...
    'beam_tilt_x_rad', tilt_x, ...
    'beam_tilt_y_rad', tilt_y);

fprintf('\nTransmitter Configuration:\n');
fprintf('  Power: %.1f W (%.1f dBm)\n', P_tx, tx_summary.P_tx_dBm);
fprintf('  Laser Linewidth: %.1f kHz\n', laser_linewidth);
fprintf('    -> Phase Noise RMS: %.3f deg\n', tx_summary.phase_noise_rms_deg);
fprintf('  Timing Jitter: %.1f ps RMS\n', timing_jitter);
fprintf('    -> Phase Error RMS: %.2f cycles (%.1f deg)\n', ...
    tx_summary.timing_jitter_phase_rms_cycles, tx_summary.timing_jitter_phase_rms_deg);
fprintf('  TX Aperture: %.1f mm radius\n', tx_aperture*1e3);
fprintf('    -> Transmission: %.2f%%\n', tx_summary.tx_aperture_transmission*100);
fprintf('    -> Aperture Loss: %.3f dB\n', tx_summary.tx_aperture_loss_dB);

%% Phase Noise Sequence Test
num_symbols = 1000;
symbol_rate = 1e9;
symbol_time = 1.0 / symbol_rate;
phase_noise_seq = beam.generate_phase_noise_sequence(num_symbols, symbol_time, laser_linewidth, 42);

fprintf('\nPhase Noise Sequence (first 5 symbols):\n');
first_5 = rad2deg(phase_noise_seq(1:5));
fprintf('  [%.4f, %.4f, %.4f, %.4f, %.4f] deg\n', first_5);

%% Visualization Logic (Equivalent to plot_beam_analysis)
% 1. Grid Setup
r_max_m = MAX_RADIUS_MM * 1e-3;
% Make grid symmetric about zero
x_vec = linspace(-r_max_m, r_max_m, GRID_SIZE);
y_vec = linspace(-r_max_m, r_max_m, GRID_SIZE);
[X, Y] = meshgrid(x_vec, y_vec);
[PHI, R_grid] = cart2pol(X, Y);

zField = 1000;
field_z0 = beam.generate_beam_field(R_grid, PHI, zField);
intensity_z0 = abs(field_z0).^2;
phase_z0 = angle(field_z0);

% 2. Radial Profile (1D)
r_range = linspace(0, r_max_m, GRID_SIZE*2);
% Use phi=0 for radial cut
intensity_profile = beam.calculate_intensity(r_range, zeros(size(r_range)), 0);

output_folder = fullfile(PLOT_DIR, sprintf('lgBeam_p%d_l%d', beam.p, beam.l));
if SAVE_FIGURE
    if ~exist(output_folder, 'dir'), mkdir(output_folder); end
end

% --- Figure 1: Radial Profile ---
fig1 = figure('Position', [100, 100, 800, 600]);
plot(r_range*1e3, intensity_profile, 'b-', 'LineWidth', 2); hold on;
xline(beam.w0*1e3, 'r--', sprintf('w_0=%.1f mm', beam.w0*1e3), 'LineWidth', 1.5);
xlabel('Radial Position r [mm]');
ylabel('Intensity [a.u.]');
title(sprintf('Lateral Intensity Profile at z = %.0f m', zField));
grid on;
if SAVE_FIGURE
    saveas(fig1, fullfile(output_folder, 'radial_profile.png'));
end

% --- Figure 2: Transverse Intensity ---
fig2 = figure('Position', [150, 150, 800, 600]);
imagesc(x_vec*1e3, y_vec*1e3, intensity_z0);
axis xy equal tight; % Important for matching Python origin='lower'
colormap('hot');
c = colorbar;
ylabel(c, 'Intensity [a.u.]');
xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Transverse Intensity at z = %.0f m', zField));
if SAVE_FIGURE
    saveas(fig2, fullfile(output_folder, 'transverse_intensity.png'));
end

% --- Figure 3: Transverse Phase ---
fig3 = figure('Position', [200, 200, 800, 600]);
imagesc(x_vec*1e3, y_vec*1e3, phase_z0);
axis xy equal tight;
% Create a cyclic colormap similar to twilight if not available
try
    colormap(twilight);
catch
    colormap(hsv); % HSV is also cyclic
end
c = colorbar;
ylabel(c, 'Phase [rad]');
clim([-pi, pi]);
set(c, 'Ticks', [-pi, 0, pi], 'TickLabels', {'-\pi', '0', '\pi'});
xlabel('x [mm]');
ylabel('y [mm]');
title(sprintf('Transverse Phase at z = %.0f m (OAM = %d)', zField, beam.l));
if SAVE_FIGURE
    saveas(fig3, fullfile(output_folder, 'transverse_phase.png'));
end

% --- Figure 4: Longitudinal Propagation ---
z_max = 3 * beam.z_R;
num_z_steps = 150;
num_r_steps = 200;

z_array = linspace(0, z_max, num_z_steps);
r_max_long = 3 * beam.beam_waist(z_max);
r_array_long = linspace(-r_max_long, r_max_long, num_r_steps);

intensity_long = zeros(num_r_steps, num_z_steps);
for i = 1:num_z_steps
    z_val = z_array(i);
    intensity_long(:, i) = beam.calculate_intensity(abs(r_array_long), zeros(size(r_array_long)), z_val);
end

w_z_array = beam.beam_waist(z_array);

fig4 = figure('Position', [250, 250, 1000, 600]);
imagesc(z_array/1000, r_array_long*1e3, intensity_long);
axis xy;
hold on;
plot(z_array/1000, w_z_array*1e3, 'c--', 'LineWidth', 2);
plot(z_array/1000, -w_z_array*1e3, 'c--', 'LineWidth', 2);
xline(beam.z_R/1000, 'g:', sprintf('z_R = %.3f km', beam.z_R/1000), 'LineWidth', 2);

colormap('hot');
c = colorbar;
ylabel(c, 'Intensity [a.u.]');
xlabel('Propagation Distance z [km]');
ylabel('Radial Position r [mm]');
title(sprintf('Longitudinal Propagation (M^2 = %d)', beam.M_squared));
if SAVE_FIGURE
    saveas(fig4, fullfile(output_folder, 'longitudinal_propagation.png'));
end

fprintf('\nAnalysis Complete. Plots saved to %s\n', output_folder);