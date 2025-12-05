% run_helico_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 633e-9;    % HeNe Red
W0_GAUSS = 1.0e-3;      % 1 mm Gaussian spot size

% Helico-Conical Parameters
% Try l=10, K=20 for a balanced spiral.
% Try l=0, K=50 for a pure "Conical" beam (circle/arc).
L_OAM = 15;             
K_PARAM = 25;           
R0_SCALE = 0.5e-3;      % 0.5 mm phase scale

% FFT Simulation Grid
% CRITICAL: High resolution is needed because the phase K*(r/r0) oscillates 
% rapidly at the edges. Low res (256/512) causes aliasing.
GRID_SIZE = 1024;       
MAX_R_IN = 2.5e-3;      % Input window radius

% Lens Setup
F_LENS = 0.4;           % 40 cm focal length

fprintf('\n--- Helico-Conical Beam Analysis ---\n');
fprintf('Parameters: l=%d, K=%.1f, r0=%.2f mm\n', L_OAM, K_PARAM, R0_SCALE*1e3);

beam = HelicoConicalBeam(L_OAM, K_PARAM, R0_SCALE, WAVELENGTH, W0_GAUSS);

%% 1. Generate Source Field (Input Plane)
vec_in = linspace(-MAX_R_IN, MAX_R_IN, GRID_SIZE);
dx_in = vec_in(2) - vec_in(1);
[X_in, Y_in] = meshgrid(vec_in, vec_in);
[PHI_in, R_in] = cart2pol(X_in, Y_in);

field_source = beam.generate_source_field(R_in, PHI_in);
int_source = abs(field_source).^2;
phase_source = angle(field_source);

%% 2. Propagate to Focus (FFT)
% The "Spiral" appears at the focal plane of a lens
field_focus = HelicoConicalBeam.propagate_to_focus(field_source);
int_focus = abs(field_focus).^2;
phase_focus = angle(field_focus);

% Calculate Physical Grid at Focus
% Fraunhofer scaling: x_focus = lambda * f * u_freq
% frequency step df = 1 / (N * dx)
% so dx_focus = lambda * f / (N * dx_in)
dx_out = (WAVELENGTH * F_LENS) / (GRID_SIZE * dx_in);
vec_out = (-(GRID_SIZE/2) : (GRID_SIZE/2 - 1)) * dx_out;

fprintf('Focal Plane Grid Resolution: %.2f um per pixel\n', dx_out*1e6);

%% Plotting

figure('Position', [100, 50, 1000, 800], 'Color', 'w');

% --- 1. Source Intensity (Gaussian) ---
subplot(2, 2, 1);
imagesc(vec_in*1e3, vec_in*1e3, int_source);
axis xy equal tight; colormap(gca, hot);
title('Source Intensity (z=0)');
xlabel('x [mm]'); ylabel('y [mm]');
% This should just look like a Gaussian blob.

% --- 2. Source Phase (The Programming) ---
subplot(2, 2, 2);
imagesc(vec_in*1e3, vec_in*1e3, phase_source);
axis xy equal tight; 
try colormap(gca, twilight); catch, colormap(gca, jet); end
colorbar;
title('Source Phase Map');
% You should see a twisted spiral phase here.

% --- 3. Focal Intensity (The Resulting Spiral) ---
subplot(2, 2, 3);
% Zoom in to central area to see the spiral clearly
zoom_factor = 0.6; 
xlim_val = vec_out(end)*1e3 * zoom_factor;

imagesc(vec_out*1e3, vec_out*1e3, int_focus);
axis xy equal tight; colormap(gca, hot);
xlim([-xlim_val, xlim_val]);
ylim([-xlim_val, xlim_val]);
title('Focal Intensity (The Spiral)');
xlabel('x_f [mm]'); ylabel('y_f [mm]');

% --- 4. 3D View of the Spiral ---
subplot(2, 2, 4);
% Downsample for surface plot performance
ds = 4; 
X_3d = vec_out(1:ds:end)*1e3;
Y_3d = vec_out(1:ds:end)*1e3;
I_3d = int_focus(1:ds:end, 1:ds:end);

% Apply crop for visualization
mask = (abs(X_3d) < xlim_val) & (abs(Y_3d) < xlim_val);
surf(X_3d, Y_3d, I_3d, 'EdgeColor', 'none');
shading interp;
colormap(gca, hot);
axis equal;
xlim([-xlim_val, xlim_val]);
ylim([-xlim_val, xlim_val]);
view(0, 90); % Top-down view
title('3D Topology of Caustic');

fprintf('Analysis Complete.\n');