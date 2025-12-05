% run_hg_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 1064e-9;   % Nd:YAG wavelength
W0 = 15e-3;             % 15 mm waist
N_MODE = 2;             % x-index (2 vertical nodes)
M_MODE = 1;             % y-index (1 horizontal node)
                        % Result: A 3x2 grid of spots

GRID_SIZE = 250;
MAX_RADIUS_MM = 60;
SAVE_FIGURE = false;
PLOT_DIR = 'plots_hg_matlab';

if ~exist(PLOT_DIR, 'dir') && SAVE_FIGURE; mkdir(PLOT_DIR); end

fprintf('\n--- Hermite-Gaussian (HG) Beam Analysis ---\n');

% Instantiate
beam = HermiteGaussianBeam(N_MODE, M_MODE, WAVELENGTH, W0);

fprintf('Beam: HG n=%d m=%d\n', beam.n, beam.m);
fprintf('  Wavelength: %.0f nm\n', beam.wavelength*1e9);
fprintf('  Waist w0: %.1f mm\n', beam.w0*1e3);
fprintf('  Rayleigh Range z0: %.2f m\n', beam.z0);
fprintf('  M^2 (x): %d, M^2 (y): %d\n', 2*beam.n+1, 2*beam.m+1);

%% Visual Setup (Cartesian Grids)
xy_max = MAX_RADIUS_MM * 1e-3;
x_vec = linspace(-xy_max, xy_max, GRID_SIZE);
y_vec = linspace(-xy_max, xy_max, GRID_SIZE);
[X, Y] = meshgrid(x_vec, y_vec);

% Generate Field at Rayleigh Range
z_field = beam.z0;
% Note: generate_beam_field accepts x,y natively
field_z = beam.generate_beam_field(X, Y, z_field, 'is_polar_input', false);
intensity_z = abs(field_z).^2;
phase_z = angle(field_z);

%% Plotting

% --- Figure 1: Transverse Intensity ---
figure('Position', [100, 100, 900, 400], 'Color', 'w');
subplot(1, 2, 1);
imagesc(x_vec*1e3, y_vec*1e3, intensity_z);
axis xy equal tight;
colormap('hot');
colorbar;
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('HG_{%d,%d} Intensity @ z = z_R', beam.n, beam.m));

% --- Figure 2: Phase (The Checkerboard) ---
subplot(1, 2, 2);
imagesc(x_vec*1e3, y_vec*1e3, phase_z);
axis xy equal tight;
try colormap(gca, twilight); catch, colormap(gca, jet); end
colorbar;
xlabel('x [mm]'); ylabel('y [mm]');
title('Phase (Notice \pi flips between lobes)');

% --- Figure 3: Diagonal 1D Cut ---
% To verify the lobes, let's cut across the diagonal
diag_coord = linspace(-xy_max, xy_max, GRID_SIZE);
% Create diagonal points
field_diag = beam.generate_beam_field(diag_coord, diag_coord, z_field);
int_diag = abs(field_diag).^2;

figure('Position', [150, 150, 600, 400], 'Color', 'w');
plot(diag_coord*1e3, int_diag, 'b-', 'LineWidth', 1.5);
grid on;
xlabel('Diagonal Position (x=y) [mm]');
ylabel('Intensity');
title(sprintf('Diagonal Cut Profile (z=%.1f m)', z_field));

