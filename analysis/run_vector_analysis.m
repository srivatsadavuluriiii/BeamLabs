% run_vector_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 1064e-9;
W0 = 2.0e-3;
GRID_SIZE = 25;         % Lower resolution for Quiver plots (arrows need space)
GRID_SIZE_HI = 200;     % High resolution for Intensity plots
MAX_R = 4.0e-3;

fprintf('\n--- Cylindrical Vector Beam Analysis ---\n');

% Create one Radial and one Azimuthal beam
beam_rad = VectorBeam('radial', WAVELENGTH, W0);
beam_azi = VectorBeam('azimuthal', WAVELENGTH, W0);

% --- Setup High Res Grids (for Intensity) ---
vec_hi = linspace(-MAX_R, MAX_R, GRID_SIZE_HI);
[X_hi, Y_hi] = meshgrid(vec_hi, vec_hi);
[PHI_hi, R_hi] = cart2pol(X_hi, Y_hi);

% --- Setup Low Res Grids (for Arrows) ---
vec_lo = linspace(-MAX_R, MAX_R, GRID_SIZE);
[X_lo, Y_lo] = meshgrid(vec_lo, vec_lo);
[PHI_lo, R_lo] = cart2pol(X_lo, Y_lo);

%% 1. Generate Fields
% Radial
[Ex_r_lo, Ey_r_lo, ~] = beam_rad.generate_vector_fields(R_lo, PHI_lo, 0);
I_rad_hi = beam_rad.calculate_intensity(R_hi, PHI_hi, 0);

% Azimuthal
[Ex_a_lo, Ey_a_lo, ~] = beam_azi.generate_vector_fields(R_lo, PHI_lo, 0);
I_azi_hi = beam_azi.calculate_intensity(R_hi, PHI_hi, 0);

%% Plotting

figure('Position', [100, 100, 1000, 500], 'Color', 'w');

% --- Radial Polarization ---
subplot(1, 2, 1);
% 1. Plot Intensity as background
imagesc(vec_hi*1e3, vec_hi*1e3, I_rad_hi);
axis xy equal tight;
colormap(gca, hot); % Black/Red background
hold on;

% 2. Plot Vector Arrows (Real part of E field at t=0)
% We normalize arrow length for visibility
Ex_real = real(Ex_r_lo);
Ey_real = real(Ey_r_lo);
norm_fac = sqrt(Ex_real.^2 + Ey_real.^2);
mask = norm_fac > max(norm_fac(:))*0.1; % Don't plot arrows in the dark center

quiver(X_lo(mask)*1e3, Y_lo(mask)*1e3, ...
       Ex_real(mask)./norm_fac(mask), Ey_real(mask)./norm_fac(mask), ...
       0.5, 'c', 'LineWidth', 1.5);

title('Radially Polarized Beam');
xlabel('x [mm]'); ylabel('y [mm]');

% --- Azimuthal Polarization ---
subplot(1, 2, 2);
imagesc(vec_hi*1e3, vec_hi*1e3, I_azi_hi);
axis xy equal tight;
colormap(gca, hot);
hold on;

Ex_real_a = real(Ex_a_lo);
Ey_real_a = real(Ey_a_lo);
norm_fac_a = sqrt(Ex_real_a.^2 + Ey_real_a.^2);
mask_a = norm_fac_a > max(norm_fac_a(:))*0.1;

quiver(X_lo(mask_a)*1e3, Y_lo(mask_a)*1e3, ...
       Ex_real_a(mask_a)./norm_fac_a(mask_a), Ey_real_a(mask_a)./norm_fac_a(mask_a), ...
       0.5, 'g', 'LineWidth', 1.5);

title('Azimuthally Polarized Beam');
xlabel('x [mm]'); ylabel('y [mm]');

sgtitle('Cylindrical Vector Beams (Arrows = Polarization Direction)');

%% 2. Through-Analyzer Simulation
% A common lab test: Place a Linear Polarizer (Analyzer) in front of the beam.
% If you pass a Radial beam through a Horizontal Polarizer, you get a "Two Lobe" pattern.

% Define Analyzer Angle (Horizontal = 0)
theta_pol = 0; 

% Project field onto polarizer axis: E_out = (E . p_hat) * p_hat
% Here we just want intensity: I_out = |Ex * cos(theta) + Ey * sin(theta)|^2

[Ex_rad_hi, Ey_rad_hi, ~] = beam_rad.generate_vector_fields(R_hi, PHI_hi, 0);
[Ex_azi_hi, Ey_azi_hi, ~] = beam_azi.generate_vector_fields(R_hi, PHI_hi, 0);

I_rad_filtered = abs(Ex_rad_hi * cos(theta_pol) + Ey_rad_hi * sin(theta_pol)).^2;
I_azi_filtered = abs(Ex_azi_hi * cos(theta_pol) + Ey_azi_hi * sin(theta_pol)).^2;

figure('Position', [150, 200, 800, 400], 'Color', 'w');
subplot(1, 2, 1);
imagesc(vec_hi*1e3, vec_hi*1e3, I_rad_filtered);
axis xy equal tight; colormap(gray);
title('Radial Beam through Horiz. Polarizer');
xlabel('x [mm]');

subplot(1, 2, 2);
imagesc(vec_hi*1e3, vec_hi*1e3, I_azi_filtered);
axis xy equal tight; colormap(gray);
title('Azimuthal Beam through Horiz. Polarizer');
xlabel('x [mm]');

fprintf('Analysis Complete.\n');