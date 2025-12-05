% run_pearcey_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 633e-9;   
W0 = 60e-6;            % Scale of the pattern (60 microns)

GRID_SIZE = 300;
MAX_X = 1.2e-3;        % View window

% Propagation Params
Z_MAX = 0.15;          % 15 cm propagation
NUM_Z = 100;

fprintf('\n--- Pearcey Beam Analysis (Cusp Catastrophe) ---\n');

% Instantiate
beam = PearceyBeam(WAVELENGTH, W0, 0);

% Grid Setup
vec = linspace(-MAX_X, MAX_X, GRID_SIZE);
dx = vec(2) - vec(1);
[X, Y] = meshgrid(vec, vec);

%% 1. Generate Source Field (The Pearcey Pattern)
fprintf('Calculating Pearcey Integral at Source...\n');
field_z0 = beam.generate_beam_field(X, Y, 0);
int_z0 = abs(field_z0).^2;
phase_z0 = angle(field_z0);

%% 2. FFT Propagation (Auto-focusing Demonstration)
z_vec = linspace(0, Z_MAX, NUM_Z);

% We will take a slice along x=0 to see the auto-focusing in the Y-Z plane
% (The Pearcey pattern is symmetric in X, but asymmetric in Y)
slice_idx = round(GRID_SIZE/2); 
int_long = zeros(GRID_SIZE, NUM_Z);

fprintf('Propagating via FFT...\n');
for i = 1:NUM_Z
    z_step = z_vec(i);
    f_prop = PearceyBeam.propagate_via_fft(field_z0, WAVELENGTH, dx, z_step);
    
    % Store the vertical cut (Y-axis)
    % The Pearcey cusp points along the Y-axis
    int_long(:, i) = abs(f_prop(:, slice_idx)).^2;
end

%% Plotting

% --- Figure 1: The Pearcey Pattern (Transverse) ---
figure('Position', [100, 100, 900, 400], 'Color', 'w');

subplot(1, 2, 1);
imagesc(vec*1e3, vec*1e3, int_z0);
axis xy equal tight; colormap(gca, hot);
title('Pearcey Intensity (z=0)');
xlabel('x [mm]'); ylabel('y [mm]');
% Note: Look for the "check mark" or "cusp" shape.

subplot(1, 2, 2);
imagesc(vec*1e3, vec*1e3, phase_z0);
axis xy equal tight; 
try colormap(gca, twilight); catch, colormap(gca, jet); end
title('Pearcey Phase');
colorbar;

% --- Figure 2: Auto-Focusing (Longitudinal) ---
figure('Position', [150, 200, 800, 500], 'Color', 'w');
imagesc(z_vec*100, vec*1e3, int_long);
axis xy;
colormap(gca, hot);
xlabel('Propagation Z [cm]');
ylabel('Y-Position [mm]');
title('Pearcey Beam Auto-Focusing (Y-Z Slice)');
colorbar;

fprintf('Analysis Complete.\n');
