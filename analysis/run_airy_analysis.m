% run_airy_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% 1. Configuration
WAVELENGTH = 532e-9;     
X0 = 50e-6;              
A_PARAM = 0.05;          

GRID_SIZE = 300;
MAX_Z_CM = 10;           % Propagation distance
PLOT_DIR = 'plots_airy_matlab';
SAVE_FIGURE = false;

if ~exist(PLOT_DIR, 'dir') && SAVE_FIGURE; mkdir(PLOT_DIR); end

fprintf('\n--- Finite-Energy Airy Beam Analysis (Optimized) ---\n');

beam = AiryBeam(WAVELENGTH, X0, A_PARAM);
summary = beam.get_beam_summary();

fprintf('  Diffraction Length z_diff: %.2f cm\n', beam.z0_airy*100);
fprintf('  Deflection at z_diff: %.1f um\n', summary.deflection_at_z0*1e6);

%% 2. Transverse Plane Calculation (z=0)
% We calculate the complex field once, then derive Intensity and Phase.
xy_max = 0.4e-3;
vec_t = linspace(-xy_max, xy_max, GRID_SIZE);
[X, Y] = meshgrid(vec_t, vec_t);

field_z0 = beam.generate_beam_field(X, Y, 0);
int_z0   = abs(field_z0).^2;
phase_z0 = angle(field_z0);

%% 3. Longitudinal Plane Calculation (Side View)
% We perform a SINGLE propagation loop to get complex field vs Z.
z_vec = linspace(0, MAX_Z_CM*1e-2, 250);
% Asymmetric x-grid to capture the upward curving trajectory
x_side_vec = linspace(-4*beam.x0, 15*beam.x0, 300); 

field_side = zeros(length(x_side_vec), length(z_vec));

fprintf('Calculating Longitudinal Propagation (Intensity + Phase)...\n');
for i = 1:length(z_vec)
    % 1D Sheet mode (is_1d=true) is faster and cleaner for side views
    field_side(:, i) = beam.generate_beam_field(x_side_vec, zeros(size(x_side_vec)), z_vec(i), 'is_1d', true);
end

int_side   = abs(field_side).^2;
phase_side = angle(field_side);
x_traj     = beam.trajectory(z_vec);

%% 4. Visualization

% --- Figure 1: Transverse Physics (The Source) ---
fig1 = figure('Position', [100, 100, 1000, 400], 'Color', 'w');

% Subplot 1: Intensity
ax1 = subplot(1, 2, 1);
imagesc(vec_t*1e3, vec_t*1e3, int_z0);
axis xy equal tight;
colormap(ax1, hot); 
colorbar;
xlabel('x [mm]'); ylabel('y [mm]');
title('Intensity (z=0): The "Corner"');

% Subplot 2: Phase
ax2 = subplot(1, 2, 2);
imagesc(vec_t*1e3, vec_t*1e3, phase_z0);
axis xy equal tight;
% Attempt cyclic colormap, fallback to jet
try colormap(ax2, twilight); catch, colormap(ax2, jet); end
clim(ax2, [-pi, pi]);
cb = colorbar; set(cb, 'Ticks', [-pi, 0, pi], 'TickLabels', {'-\pi', '0', '\pi'});
xlabel('x [mm]'); ylabel('y [mm]');
title('Phase (z=0): The "Checkerboard"');


% --- Figure 2: Longitudinal Physics (The Acceleration) ---
fig2 = figure('Position', [150, 150, 1000, 600], 'Color', 'w');

% Subplot 1: Intensity Trajectory
ax3 = subplot(2, 1, 1);
imagesc(z_vec*100, x_side_vec*1e3, int_side);
axis xy;
colormap(ax3, flipud(gray)); % Inverted gray for contrast
hold on;
plot(z_vec*100, x_traj*1e3, 'r--', 'LineWidth', 2);
ylabel('x [mm]');
title(sprintf('Intensity Acceleration (a=%.2f)', beam.a));
legend('Theory (x \propto z^2)', 'Location', 'best');

% Subplot 2: Phase Wavefronts
ax4 = subplot(2, 1, 2);
imagesc(z_vec*100, x_side_vec*1e3, phase_side);
axis xy;
try colormap(ax4, twilight); catch, colormap(ax4, hsv); end
hold on;
plot(z_vec*100, x_traj*1e3, 'k--', 'LineWidth', 1.5); % Black line for visibility on color
xlabel('Propagation Z [cm]');
ylabel('x [mm]');
title('Wavefront Bending (Phase)');

fprintf('Analysis Complete.\n');