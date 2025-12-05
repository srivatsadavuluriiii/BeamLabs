% run_mathieu_analysis.m
clear; clc; close all;

% Add parent directory to path to access beam classes
script_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(script_dir);
addpath(parent_dir);

%% Configuration
WAVELENGTH = 633e-9;   % HeNe Red
W0_APERTURE = 50e-6;   % 50 micron spot size scale

% Mathieu Parameters
ORDER_M = 2;

% Q parameter controls ellipticity.
Q_LIST = [0, 5, 20];
PARITY = 'even'; % 'even' (ce) or 'odd' (se)

GRID_SIZE = 250;
MAX_R = 80e-6; 

fprintf('\n--- Mathieu-Gaussian Beam Analysis ---\n');
fprintf('Order m=%d, Parity=%s\n', ORDER_M, PARITY);

vec = linspace(-MAX_R, MAX_R, GRID_SIZE);
[X, Y] = meshgrid(vec, vec);

%% 1. Q-Parameter Transition
figure('Position', [50, 100, 1200, 350], 'Color', 'w');

for i = 1:length(Q_LIST)
    q_val = Q_LIST(i);
    
    % Create Beam
    beam = MathieuBeam(ORDER_M, q_val, PARITY, WAVELENGTH, W0_APERTURE);
    
    % Generate Field
    field = beam.generate_beam_field(X, Y, 0);
    intensity = abs(field).^2;
    
    % Plot
    subplot(1, 3, i);
    imagesc(vec*1e6, vec*1e6, intensity);
    axis xy equal tight;
    colormap(hot);
    title(sprintf('q=%.1f (m=%d)', q_val, ORDER_M));
    xlabel('x [\mum]'); ylabel('y [\mum]');
end
sgtitle('Mathieu Beam Shape Evolution');

%% 2. Detailed Analysis (q = 25)
Q_TARGET = 25.0;
beam_main = MathieuBeam(ORDER_M, Q_TARGET, PARITY, WAVELENGTH, W0_APERTURE);

fprintf('\nCalculated for q=%.1f:\n', Q_TARGET);
fprintf('  Transverse k_t: %.2e rad/m\n', beam_main.kt);
fprintf('  Semifocal distance f: %.2f um\n', beam_main.f_focal*1e6);

f_main = beam_main.generate_beam_field(X, Y, 0);
int_main = abs(f_main).^2;
phase_main = angle(f_main);

figure('Position', [100, 200, 800, 400], 'Color', 'w');
subplot(1, 2, 1);
imagesc(vec*1e6, vec*1e6, int_main);
axis xy equal tight; colormap(gca, hot);
title(sprintf('Intensity (q=%.1f)', Q_TARGET));

subplot(1, 2, 2);
imagesc(vec*1e6, vec*1e6, phase_main);
axis xy equal tight; 
try colormap(gca, twilight); catch, colormap(gca, jet); end
title('Phase');
colorbar;

fprintf('Analysis Complete.\n');