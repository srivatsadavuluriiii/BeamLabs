% run_flatTop_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 1064e-9;   % YAG Laser
W0 = 1.5e-3;            % 1.5 mm radius

% Comparison List: 0 = Gaussian, High = Sharp Top-Hat
ORDERS_TO_COMPARE = [0, 1, 2, 3, 4, 5]; 
TARGET_ORDER_FOR_PROP = 5; % Which one to simulate propagating

GRID_SIZE_1D = 300;
GRID_SIZE_2D = 200;     % Lower res for 2D maps to be fast
MAX_R = 3.0e-3;         % View window

fprintf('\n--- Flat-Top (Flattened Gaussian) Beam Analysis ---\n');

%% 1. Multi-Order Comparison (Source Plane z=0)
r_vec = linspace(0, MAX_R, GRID_SIZE_1D);
phi_dummy = zeros(size(r_vec));

% Prepare 2D Grid
vec_2d = linspace(-MAX_R, MAX_R, GRID_SIZE_2D);
[X, Y] = meshgrid(vec_2d, vec_2d);
[PHI_2D, R_2D] = cart2pol(X, Y);

% Storage for plotting
profiles_1d = zeros(length(ORDERS_TO_COMPARE), GRID_SIZE_1D);
maps_2d = cell(1, length(ORDERS_TO_COMPARE));

fprintf('Calculating profiles for N = %s...\n', num2str(ORDERS_TO_COMPARE));

for i = 1:length(ORDERS_TO_COMPARE)
    n_val = ORDERS_TO_COMPARE(i);
    beam_tmp = FlatTopBeam(n_val, W0, WAVELENGTH);
    
    % 1D Cut
    I_1d = beam_tmp.calculate_intensity(r_vec, phi_dummy, 0);
    profiles_1d(i, :) = I_1d / max(I_1d); % Normalize
    
    % 2D Map
    field_2d = beam_tmp.generate_beam_field(R_2D, PHI_2D, 0);
    I_2d = abs(field_2d).^2;
    maps_2d{i} = I_2d / max(I_2d(:)); % Normalize
end

%% 2. Longitudinal Evolution (Target Order Only)
% Simulating the "Hot Edge" diffraction for the sharpest beam
beam_prop = FlatTopBeam(TARGET_ORDER_FOR_PROP, W0, WAVELENGTH);
z_max = 0.5 * beam_prop.z0; 
z_vec = linspace(0, z_max, 250);
r_prop = linspace(-2*W0, 2*W0, 200);

int_long = zeros(length(r_prop), length(z_vec));

fprintf('Calculating Propagation for N=%d...\n', TARGET_ORDER_FOR_PROP);
for i = 1:length(z_vec)
    f = beam_prop.generate_beam_field(abs(r_prop), zeros(size(r_prop)), z_vec(i));
    int_long(:, i) = abs(f).^2;
end

%% Plotting

% --- Figure 1: 1D Profile Overlays ---
figure('Position', [50, 50, 600, 400], 'Color', 'w');
hold on;
colors = lines(length(ORDERS_TO_COMPARE));
legend_str = cell(1, length(ORDERS_TO_COMPARE));

for i = 1:length(ORDERS_TO_COMPARE)
    plot(r_vec*1e3, profiles_1d(i, :), 'LineWidth', 2, 'Color', colors(i,:));
    legend_str{i} = sprintf('N = %d', ORDERS_TO_COMPARE(i));
end

xline(W0*1e3, 'k:', 'Waist w0');
legend(legend_str, 'Location', 'SouthWest');
xlabel('Radial Position r [mm]');
ylabel('Normalized Intensity');
title('Transition from Gaussian to Flat-Top');
grid on;

% --- Figure 2: 2D Intensity Maps (Side-by-Side) ---
figure('Position', [100, 100, 1200, 300], 'Color', 'w');
for i = 1:length(ORDERS_TO_COMPARE)
    subplot(1, length(ORDERS_TO_COMPARE), i);
    imagesc(vec_2d*1e3, vec_2d*1e3, maps_2d{i});
    axis xy equal tight;
    colormap(hot);
    title(sprintf('N = %d', ORDERS_TO_COMPARE(i)));
    xlabel('x [mm]'); 
    if i==1, ylabel('y [mm]'); end
end
sgtitle('Transverse Beam Shapes (z=0)');

% --- Figure 3: 3D Surface Plot (High Order) ---
% This helps visualize the "Table Top" flatness
figure('Position', [150, 150, 600, 500], 'Color', 'w');
% Use the last map (highest N)
surf(X*1e3, Y*1e3, maps_2d{end}, 'EdgeColor', 'none');
shading interp;
colormap(jet);
lightangle(-45, 30);
lighting gouraud;
material shiny;
axis equal;
title(sprintf('3D Intensity Topology (N=%d)', ORDERS_TO_COMPARE(end)));
xlabel('x [mm]'); ylabel('y [mm]'); zlabel('Intensity');
view(30, 45);

% --- Figure 4: Propagation Dynamics ---
figure('Position', [200, 200, 900, 500], 'Color', 'w');
imagesc(z_vec, r_prop*1e3, int_long);
axis xy;
colormap(hot);
colorbar;
xlabel('Propagation Distance z [m]');
ylabel('Radial Position r [mm]');
title(sprintf('Diffraction of Flat-Top (N=%d)', TARGET_ORDER_FOR_PROP));

hold on;
yline(W0*1e3, 'c--', 'LineWidth', 1);
yline(-W0*1e3, 'c--', 'LineWidth', 1);
text(z_max*0.05, W0*1e3*1.1, 'Beam Edge', 'Color', 'c');

fprintf('Analysis Complete.\n');