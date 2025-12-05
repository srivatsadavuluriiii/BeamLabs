% run_ince_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 1064e-9;
W0 = 2.0e-3;

% Ince Parameters
P_ORDER = 5;
M_DEGREE = 3; 
% Epsilon controls the "Ellipticity"
% 0.5 = Near Circular (LG-like)
% 2.0 = Elliptical (True Ince)
% 8.0 = Rectangular (HG-like)
EPSILON_LIST = [0.5, 2.0, 8.0]; 
PARITY = 'even';

GRID_SIZE = 250;
MAX_R = 6.0e-3; % View window radius

fprintf('\n--- Ince-Gaussian Beam Analysis ---\n');
fprintf('Order p=%d, Degree m=%d, Parity=%s\n', P_ORDER, M_DEGREE, PARITY);

vec = linspace(-MAX_R, MAX_R, GRID_SIZE);
[X, Y] = meshgrid(vec, vec);

%% 1. Epsilon Comparison (Shape Evolution)
fig1 = figure('Position', [50, 100, 1200, 350], 'Color', 'w');

for i = 1:length(EPSILON_LIST)
    eps_val = EPSILON_LIST(i);
    
    % Create Beam
    beam = InceGaussianBeam(P_ORDER, M_DEGREE, eps_val, PARITY, WAVELENGTH, W0);
    
    % Generate
    field = beam.generate_beam_field(X, Y, 0);
    intensity = abs(field).^2;
    
    % Plot
    subplot(1, 3, i);
    imagesc(vec*1e3, vec*1e3, intensity);
    axis xy equal tight;
    colormap(hot);
    title(sprintf('\\epsilon = %.1f (Transition)', eps_val));
    xlabel('x [mm]'); ylabel('y [mm]');
end
sgtitle('Transition from LG-like (Left) to HG-like (Right)');

%% 2. Detailed Analysis for Target Epsilon
EPS_TARGET = 2.0; % The "True" Ince regime
beam_main = InceGaussianBeam(P_ORDER, M_DEGREE, EPS_TARGET, PARITY, WAVELENGTH, W0);

fprintf('\nDetailed Analysis for epsilon = %.1f:\n', EPS_TARGET);
fprintf('  Rayleigh Range z0: %.2f m\n', beam_main.z0);

% --- Generate Fields ---
field_main = beam_main.generate_beam_field(X, Y, 0);
int_main = abs(field_main).^2;
phase_main = angle(field_main);

% --- 1D Cuts (Major vs Minor Axis) ---
% IG beams are not radially symmetric, so x and y profiles differ
x_cut_int = beam_main.calculate_intensity(vec, zeros(size(vec)), 0);
y_cut_int = beam_main.calculate_intensity(zeros(size(vec)), vec, 0);

% --- Longitudinal Propagation (X-Z Slice) ---
z_max = 2.0 * beam_main.z0;
num_z = 150;
z_vec = linspace(0, z_max, num_z);
x_prop_vec = linspace(-2.5*W0, 2.5*W0, 200); 

% Pre-allocate
int_long = zeros(length(x_prop_vec), length(z_vec));

fprintf('Calculating Longitudinal Propagation...\n');
for i = 1:length(z_vec)
    % Slice along y=0 (Major Axis view)
    f_step = beam_main.generate_beam_field(x_prop_vec, zeros(size(x_prop_vec)), z_vec(i));
    int_long(:, i) = abs(f_step).^2;
end

%% Plotting Detailed Views

% --- Figure 2: Transverse Physics (Intensity & Phase) ---
figure('Position', [100, 200, 900, 400], 'Color', 'w');

% Intensity Map
subplot(1, 2, 1);
imagesc(vec*1e3, vec*1e3, int_main);
axis xy equal tight; 
colormap(gca, hot);
xlabel('x [mm]'); ylabel('y [mm]');
title(sprintf('Transverse Intensity (\\epsilon=%.1f)', EPS_TARGET));
colorbar;

% Phase Map
subplot(1, 2, 2);
imagesc(vec*1e3, vec*1e3, phase_main);
axis xy equal tight; 
try colormap(gca, twilight); catch, colormap(gca, jet); end
xlabel('x [mm]'); ylabel('y [mm]');
title('Transverse Phase Structure');
cb = colorbar; set(cb, 'Ticks', [-pi, 0, pi], 'TickLabels', {'-\pi', '0', '\pi'});

% --- Figure 3: 1D Axial Profiles (Anisotropy) ---
figure('Position', [150, 300, 600, 400], 'Color', 'w');
plot(vec*1e3, x_cut_int, 'b-', 'LineWidth', 2); hold on;
plot(vec*1e3, y_cut_int, 'r--', 'LineWidth', 1.5);
grid on;
xlabel('Position [mm]');
ylabel('Intensity [a.u.]');
title('Axial Cross-Sections (Anisotropy Check)');
legend('Major Axis (x)', 'Minor Axis (y)');
xlim([-MAX_R*1e3, MAX_R*1e3]);

% --- Figure 4: Longitudinal Propagation (Side View) ---
figure('Position', [200, 400, 800, 500], 'Color', 'w');
imagesc(z_vec, x_prop_vec*1e3, int_long);
axis xy;
colormap(gca, hot);
xlabel('Propagation Distance z [m]');
ylabel('Transverse Position x [mm]');
title(sprintf('Longitudinal Propagation (X-Z Plane, \\epsilon=%.1f)', EPS_TARGET));
colorbar;

% Overlay beam expansion lines
w_z_prop = beam_main.w0 * sqrt(1 + (z_vec/beam_main.z0).^2);
hold on;
plot(z_vec, w_z_prop*1e3, 'c--', 'LineWidth', 1.5);
plot(z_vec, -w_z_prop*1e3, 'c--', 'LineWidth', 1.5);
xline(beam_main.z0, 'g:', 'Rayleigh Range z_R');
legend('Intensity', 'Gaussian Width w(z)', 'Rayleigh Range');

fprintf('Analysis Complete.\n');