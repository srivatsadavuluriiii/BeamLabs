% run_pov_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 1550e-9;
R_RING = 1.5e-3;        % Fixed Ring Radius: 1.5 mm
W0_THICKNESS = 0.3e-3;  % Ring Thickness: 0.3 mm

% We will compare two modes
L1 = 2;
L2 = 10; % Much higher OAM

GRID_SIZE = 300;
MAX_R = 3.0e-3;         % View window
PLOT_DIR = 'plots_pov_matlab';
if ~exist(PLOT_DIR, 'dir'); mkdir(PLOT_DIR); end

fprintf('\n--- Perfect Optical Vortex (POV) Analysis ---\n');
fprintf('Target Radius: %.2f mm (Independent of l)\n', R_RING*1e3);

% Instantiate two beams
beam1 = PerfectVortexBeam(L1, R_RING, W0_THICKNESS, WAVELENGTH);
beam2 = PerfectVortexBeam(L2, R_RING, W0_THICKNESS, WAVELENGTH);

%% 1. Transverse Profile Calculation
vec = linspace(-MAX_R, MAX_R, GRID_SIZE);
[X, Y] = meshgrid(vec, vec);
[PHI, R] = cart2pol(X, Y);

% Generate Fields at z=0
f1 = beam1.generate_beam_field(R, PHI, 0);
f2 = beam2.generate_beam_field(R, PHI, 0);

i1 = abs(f1).^2;
i2 = abs(f2).^2;

% Normalize for comparison (since power density varies with l)
i1 = i1 / max(i1(:));
i2 = i2 / max(i2(:));

% Get Phases
p1 = angle(f1);
p2 = angle(f2);

%% 2. Radial Cut Comparison (The Proof)
% Extract a slice along x-axis (phi=0)
r_cut = linspace(0, MAX_R, GRID_SIZE);
phi_cut = zeros(size(r_cut));

int_cut1 = beam1.calculate_intensity(r_cut, phi_cut, 0);
int_cut2 = beam2.calculate_intensity(r_cut, phi_cut, 0);

% Normalize cuts
int_cut1 = int_cut1 / max(int_cut1);
int_cut2 = int_cut2 / max(int_cut2);

%% Plotting

% --- Figure 1: Visual Comparison ---
figure('Position', [100, 100, 1000, 500], 'Color', 'w');

% Mode 1 Intensity
subplot(2, 2, 1);
imagesc(vec*1e3, vec*1e3, i1);
axis xy equal tight; colormap(gca, hot);
title(sprintf('POV Intensity (l=%d)', L1));
xlabel('x [mm]'); ylabel('y [mm]');

% Mode 1 Phase
subplot(2, 2, 2);
imagesc(vec*1e3, vec*1e3, p1);
axis xy equal tight; 
try colormap(gca, twilight); catch, colormap(gca, jet); end
title(sprintf('POV Phase (l=%d)', L1));
colorbar;

% Mode 2 Intensity (Should look SAME size)
subplot(2, 2, 3);
imagesc(vec*1e3, vec*1e3, i2);
axis xy equal tight; colormap(gca, hot);
title(sprintf('POV Intensity (l=%d)', L2));
xlabel('x [mm]'); ylabel('y [mm]');

% Mode 2 Phase (More spokes)
subplot(2, 2, 4);
imagesc(vec*1e3, vec*1e3, p2);
axis xy equal tight; 
try colormap(gca, twilight); catch, colormap(gca, jet); end
title(sprintf('POV Phase (l=%d)', L2));
colorbar;

% --- Figure 2: The Proof of "Perfectness" ---
figure('Position', [150, 150, 600, 400], 'Color', 'w');
plot(r_cut*1e3, int_cut1, 'b-', 'LineWidth', 2); hold on;
plot(r_cut*1e3, int_cut2, 'r--', 'LineWidth', 2);

xline(R_RING*1e3, 'k:', 'Target Radius R_r');
legend(sprintf('l = %d', L1), sprintf('l = %d', L2), 'Target Radius');

xlabel('Radial Position [mm]');
ylabel('Normalized Intensity');
title('Radial Profiles Overlay (Proof of Radius Independence)');
grid on;
xlim([0, MAX_R*1e3]);

fprintf('Analysis Complete.\n');