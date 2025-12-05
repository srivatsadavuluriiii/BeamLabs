% run_frozen_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 532e-9;    % Green
W0_APERTURE = 3.0e-3;   % 3 mm Aperture 

% Frozen Wave Parameters
L_PATTERN = 0.10;       % 10 cm Total Pattern Length
N_MODES = 30;           % Resolution of the pattern (Total Bessel beams = 2*N + 1)
Q_FACTOR = 0.9998 * (2*pi/WAVELENGTH); 

% Target: "Dash-Gap-Dash"
% z in meters: ON at [2-4]cm and [6-8]cm
target_profile = @(z) double( (z > 0.02 & z < 0.04) | (z > 0.06 & z < 0.08) );

fprintf('\n--- Frozen Wave (Longitudinal Shaping) Analysis ---\n');
beam = FrozenWaveBeam(WAVELENGTH, W0_APERTURE, L_PATTERN, Q_FACTOR, N_MODES, target_profile);

%% 1. Longitudinal Calculation (Side View X-Z Plane)
z_vec = linspace(0, L_PATTERN, 400);
x_vec = linspace(-100e-6, 100e-6, 200); 
int_long = zeros(length(x_vec), length(z_vec));

fprintf('Calculating X-Z Longitudinal Field...\n');
for i = 1:length(z_vec)
    f = beam.generate_beam_field(abs(x_vec), zeros(size(x_vec)), z_vec(i));
    int_long(:, i) = abs(f).^2;
end
int_long = int_long / max(int_long(:));

%% 2. Transverse Calculation (Slices at specific Z)
% Slice 1: Inside the first "Dash" (z = 3 cm)
z_slice_ON = 0.03; 
% Slice 2: Inside the "Gap" (z = 5 cm)
z_slice_OFF = 0.05;

xy_vec = linspace(-80e-6, 80e-6, 200);
[X, Y] = meshgrid(xy_vec, xy_vec);
[PHI, R] = cart2pol(X, Y);

fprintf('Calculating Transverse Slices...\n');
field_on  = beam.generate_beam_field(R, PHI, z_slice_ON);
field_off = beam.generate_beam_field(R, PHI, z_slice_OFF);

I_on = abs(field_on).^2;   I_on = I_on / max(I_on(:));
I_off = abs(field_off).^2; % Normalize to the ON peak to show relative darkness

%% 3. 3D Volume Generation (For Isosurface)
fprintf('Calculating 3D Volume for Rendering...\n');

% Define 3D coordinates
z_3d = linspace(0, L_PATTERN, 100);
xy_3d = linspace(-60e-6, 60e-6, 60);

% CRITICAL FIX: Create full 3D grids to match the volume data structure
[X3D, Y3D, Z3D] = meshgrid(xy_3d, xy_3d, z_3d);
[PHI3D, R3D] = cart2pol(X3D, Y3D);

vol_3d = zeros(length(xy_3d), length(xy_3d), length(z_3d));

for k = 1:length(z_3d)
    % We extract the 2D slice of coordinates for the current z
    r_slice = R3D(:, :, k);
    p_slice = PHI3D(:, :, k);
    z_val = z_3d(k);
    
    f_step = beam.generate_beam_field(r_slice, p_slice, z_val);
    vol_3d(:, :, k) = abs(f_step).^2;
end
vol_3d = vol_3d / max(vol_3d(:));

%% Plotting

% --- Figure 1: Design vs Reality (Longitudinal) ---
figure('Position', [50, 50, 800, 600], 'Color', 'w');
subplot(2, 1, 1);
target_vals = target_profile(z_vec);
center_idx = round(length(x_vec)/2);
actual_vals = int_long(center_idx, :);

plot(z_vec*100, target_vals, 'k--', 'LineWidth', 1.5); hold on;
plot(z_vec*100, actual_vals, 'r-', 'LineWidth', 1.5);
legend('Target F(z)', 'Simulated');
xlabel('z [cm]'); ylabel('Intensity');
title('Longitudinal On-Axis Profile');
grid on; ylim([-0.1, 1.1]);

subplot(2, 1, 2);
imagesc(z_vec*100, x_vec*1e6, int_long);
axis xy; colormap(hot);
xlabel('z [cm]'); ylabel('r [\mum]');
title('Side View (X-Z Plane)');

% --- Figure 2: The "Recipe" (Spectral Weights) ---
figure('Position', [100, 100, 600, 400], 'Color', 'w');
n_indices = -N_MODES:N_MODES;
stem(n_indices, abs(beam.A_coeffs), 'filled', 'BaseValue', 0, 'MarkerFaceColor', 'b');
xlabel('Mode Index n (Bessel Order)');
ylabel('Magnitude |A_n|');
title('Spectral Weights (The Fourier Transform of the Pattern)');
grid on;

% --- Figure 3: Transverse Slices (ON vs OFF) ---
figure('Position', [150, 150, 800, 400], 'Color', 'w');

subplot(1, 2, 1);
imagesc(xy_vec*1e6, xy_vec*1e6, I_on);
axis xy equal tight; colormap(gca, hot); colorbar;
title(sprintf('Slice @ z=%.1f cm (ON)', z_slice_ON*100));
xlabel('x [\mum]'); ylabel('y [\mum]');

subplot(1, 2, 2);
imagesc(xy_vec*1e6, xy_vec*1e6, I_off);
axis xy equal tight; colormap(gca, hot); colorbar;
title(sprintf('Slice @ z=%.1f cm (OFF)', z_slice_OFF*100));
xlabel('x [\mum]');

% --- Figure 4: 3D Light Bullets ---
figure('Position', [200, 200, 700, 600], 'Color', 'k');

% Use the consistent 3D Grids (X3D, Y3D, Z3D)
p = patch(isosurface(X3D*1e6, Y3D*1e6, Z3D*100, vol_3d, 0.4));

set(p, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.8);
isonormals(X3D*1e6, Y3D*1e6, Z3D*100, vol_3d, p);

view(35, 30); 
axis tight; camlight; lighting gouraud;
xlabel('x [\mum]'); ylabel('y [\mum]'); zlabel('z [cm]');
title('3D Iso-surface (Light Bullets)', 'Color', 'w');
set(gca, 'XColor', 'w', 'YColor', 'w', 'ZColor', 'w');
grid on;

fprintf('Analysis Complete.\n');