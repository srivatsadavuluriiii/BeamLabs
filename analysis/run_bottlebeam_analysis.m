% run_bottlebeam_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 532e-9;    % Green trapping laser
W0 = 2.0e-3;            % 2 mm aperture

% Bottle Parameters
% beta determines the radial size
% delta_beta determines the longitudinal length
BETA_CENTER = 30000;    % High beta = small radial bottle
DELTA_BETA  = 2000;     % Difference between the two beams

GRID_SIZE = 200;
MAX_R_UM = 150;         % View window (microns)
MAX_Z_MM = 10;          % Propagation window (mm)

fprintf('\n--- Optical Bottle Beam Analysis ---\n');

beam = BottleBeam(WAVELENGTH, BETA_CENTER, DELTA_BETA, W0);

fprintf('Parameters:\n');
fprintf('  Radial Wavevector (beta): %.0f rad/m\n', beam.beta_r);
fprintf('  Bottle Length Parameter (L): %.2f mm\n', beam.delta_z*1e3);
fprintf('  Max Range (Gaussian limit): %.2f mm\n', beam.z_max*1e3);

%% 1. Transverse View (The "Mouth" of the Bottle) at z=0
% Ideally, this should be DARK in the center, surrounded by a bright ring.
vec_r = linspace(-MAX_R_UM*1e-6, MAX_R_UM*1e-6, GRID_SIZE);
[X, Y] = meshgrid(vec_r, vec_r);
[PHI, R] = cart2pol(X, Y);

field_z0 = beam.generate_beam_field(R, PHI, 0);
int_z0 = abs(field_z0).^2;

%% 2. Longitudinal View (The "Bottle" Shape)
% We slice through X-Z to see the 3D void.
z_vec = linspace(-MAX_Z_MM*1e-3, MAX_Z_MM*1e-3, 300);
x_vec_long = linspace(-MAX_R_UM*1e-6, MAX_R_UM*1e-6, 200);

int_long = zeros(length(x_vec_long), length(z_vec));

fprintf('Calculating Bottle Geometry...\n');
for i = 1:length(z_vec)
    % 1D slice along x
    f = beam.generate_beam_field(abs(x_vec_long), zeros(size(x_vec_long)), z_vec(i));
    int_long(:, i) = abs(f).^2;
end

%% Plotting

% --- Figure 1: Transverse Cross-Section at z=0 ---
figure('Position', [100, 100, 500, 400], 'Color', 'w');
imagesc(vec_r*1e6, vec_r*1e6, int_z0);
axis xy equal tight;
colormap(hot);
colorbar;
xlabel('x [\mum]'); ylabel('y [\mum]');
title('Transverse Cross-Section (z=0)');
% Note: You should see a dark hole in the center.

% --- Figure 2: Longitudinal "Bottle" Structure ---
figure('Position', [200, 150, 900, 500], 'Color', 'w');
imagesc(z_vec*1e3, x_vec_long*1e6, int_long);
axis xy;
colormap(hot);
colorbar;
xlabel('Propagation z [mm]');
ylabel('Radial position r [\mum]');
title('Optical Bottle Structure (Side View)');

hold on;
% Verify the "walls"
yline(0, 'w:', 'Axis');
xline(0, 'w:', 'Focus');

% Add contour to highlight the void
contour(z_vec*1e3, x_vec_long*1e6, int_long, 10, 'c', 'LineWidth', 0.5);

% Mark theoretical bottle length
L_theory = beam.delta_z;
line([-L_theory/2, L_theory/2]*1e3, [0, 0], 'Color', 'g', 'LineWidth', 2, 'LineStyle', '--');
text(0, -MAX_R_UM*0.8, 'Theoretical Length', 'Color', 'g', 'HorizontalAlignment', 'center');

fprintf('Analysis Complete.\n');