# BeamLabs

A comprehensive MATLAB library for modeling and analyzing structured light beams. BeamLabs provides implementations of 14 different beam types with consistent interfaces, advanced propagation physics, and detailed analysis capabilities.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Beam Classes](#beam-classes)
- [Analysis Scripts](#analysis-scripts)
- [Common Interface](#common-interface)
- [Advanced Features](#advanced-features)
- [Examples](#examples)
- [References](#references)

## Overview

BeamLabs is designed for researchers and engineers working with structured light in applications such as:
- Free-space optical communications
- Optical trapping and manipulation
- Quantum optics
- Beam shaping and propagation analysis
- Laser mode characterization

The library implements both fundamental beams (Gaussian, Laguerre-Gaussian, Hermite-Gaussian) and specialized beams (Airy, Bessel-Gauss, Perfect Vortex, Vector beams, etc.) with rigorous physical modeling.

## Features

- **14 Beam Types**: From fundamental Gaussian to exotic Pearcey and Mathieu beams
- **Consistent Interface**: All beam classes share common methods and parameter conventions
- **Advanced Physics**: Gouy phase, wavefront curvature, beam waist evolution, Rayleigh range
- **Noise Modeling**: Phase noise and timing jitter for realistic simulations
- **Propagation Analysis**: Longitudinal propagation visualization and parameter tracking
- **Polarization Support**: Vector beams with spatially varying polarization
- **Numerical Solvers**: Custom implementations for Ince, Mathieu, and Pearcey functions
- **Comprehensive Analysis**: Pre-built analysis scripts for each beam type

## Installation

1. Clone or download this repository
2. Add the BeamLabs directory to your MATLAB path:
   ```matlab
   addpath('/path/to/BeamLabs');
   ```
3. Ensure you have MATLAB R2018b or later (for classdef support)

**Required Toolboxes**: None (all implementations use base MATLAB)

## Quick Start

### Basic Example: Gaussian Beam

```matlab
% Create a Gaussian beam
wavelength = 1550e-9;  % 1550 nm (telecom C-band)
w0 = 2e-3;            % 2 mm beam waist
beam = GaussianBeam(wavelength, w0);

% Generate field at z=0
r = linspace(0, 5e-3, 200);
phi = zeros(size(r));
z = 0;
field = beam.generate_beam_field(r, phi, z);

% Calculate intensity
intensity = abs(field).^2;

% Plot
plot(r*1e3, intensity);
xlabel('Radius [mm]');
ylabel('Intensity [a.u.]');
title('Gaussian Beam Profile');
```

### Example: Laguerre-Gaussian with OAM

```matlab
% Create LG beam with OAM
wavelength = 1550e-9;
w0 = 25e-3;
p = 0;  % Radial index
l = 2;  % Azimuthal index (OAM)
beam = LaguerreGaussianBeam(p, l, wavelength, w0);

% Generate 2D field
grid_size = 200;
r_max = 75e-3;
x = linspace(-r_max, r_max, grid_size);
y = linspace(-r_max, r_max, grid_size);
[X, Y] = meshgrid(x, y);
[PHI, R] = cart2pol(X, Y);

field = beam.generate_beam_field(R, PHI, 0);
intensity = abs(field).^2;

% Visualize
imagesc(x*1e3, y*1e3, intensity);
axis xy equal tight;
colormap('hot');
colorbar;
```

## Beam Classes

### 1. GaussianBeam
**Fundamental TEM00 paraxial Gaussian beam**

```matlab
beam = GaussianBeam(wavelength, w0)
```

**Parameters:**
- `wavelength`: Wavelength in meters
- `w0`: Beam waist radius at z=0 (1/e² intensity)

**Key Methods:**
- `beam_waist(z)`: Calculate beam waist at distance z
- `radius_of_curvature(z)`: Wavefront curvature
- `gouy_phase(z)`: Gouy phase shift
- `calculate_coupling_loss(r_aper)`: Aperture transmission loss

**Use Cases:** Baseline for optical communications, laser physics

---

### 2. LaguerreGaussianBeam
**LG modes with orbital angular momentum (OAM)**

```matlab
beam = LaguerreGaussianBeam(p, l, wavelength, w0)
```

**Parameters:**
- `p`: Radial mode index (≥0)
- `l`: Azimuthal mode index (OAM topological charge)
- `wavelength`: Wavelength in meters
- `w0`: Beam waist radius

**Key Features:**
- OAM phase structure: `exp(-i*l*phi)`
- M² beam quality factor calculation
- Phase noise and timing jitter modeling
- Aperture loss calculation
- Overlap integrals with other beams

**Use Cases:** Free-space optical communications, quantum optics, optical tweezers

---

### 3. HermiteGaussianBeam
**HG modes in Cartesian symmetry**

```matlab
beam = HermiteGaussianBeam(n, m, wavelength, w0)
```

**Parameters:**
- `n`: Mode index along x-axis (≥0)
- `m`: Mode index along y-axis (≥0)
- `wavelength`: Wavelength in meters
- `w0`: Beam waist radius

**Key Features:**
- Cartesian mode patterns
- Hermite polynomial evaluation
- Separate M² for x and y axes

**Use Cases:** Laser mode analysis, resonator characterization

---

### 4. BesselGaussianBeam
**Bessel-Gauss beams with diffraction-free properties**

```matlab
beam = BesselGaussianBeam(l, beta, wavelength, w0)
```

**Parameters:**
- `l`: Topological charge (Bessel order)
- `beta`: Transverse wavenumber (rad/m) - controls ring spacing
- `wavelength`: Wavelength in meters
- `w0`: Gaussian envelope waist

**Key Features:**
- Diffraction-free range calculation (`z_max`)
- Geometric cone angle
- Propagation summary (Bessel zone vs. far field)
- Phase noise support

**Use Cases:** Non-diffracting beams, optical trapping, long-range propagation

---

### 5. AiryBeam
**Finite-energy Airy beams with parabolic acceleration**

```matlab
beam = AiryBeam(wavelength, x0, a)
```

**Parameters:**
- `wavelength`: Wavelength in meters
- `x0`: Transverse scale parameter (main lobe width)
- `a`: Apodization parameter (0 < a << 1, controls decay)

**Key Features:**
- Parabolic trajectory: `x ∝ z²`
- 1D sheet mode option (`is_1d` parameter)
- Trajectory calculation method
- Characteristic diffraction length

**Use Cases:** Self-accelerating beams, curved light paths, optical manipulation

---

### 6. BottleBeam
**Optical bottle (dark focus surrounded by bright regions)**

```matlab
beam = BottleBeam(wavelength, beta_center, delta_beta, w0)
```

**Parameters:**
- `wavelength`: Wavelength in meters
- `beta_center`: Average radial wavenumber
- `delta_beta`: Difference in radial wavenumbers (controls bottle length)
- `w0`: Gaussian apodization waist

**Key Features:**
- Superposition of two Bessel beams
- Destructive interference at center
- Configurable bottle length

**Use Cases:** Optical trapping, particle manipulation, dark beam generation

---

### 7. FlatTopBeam
**Flattened Gaussian (top-hat) beam**

```matlab
beam = FlatTopBeam(N_order, w0, wavelength)
```

**Parameters:**
- `N_order`: Flatness order (0 = Gaussian, 20+ = sharp top-hat)
- `w0`: Waist parameter
- `wavelength`: Wavelength in meters

**Key Features:**
- Gori's Flattened Gaussian expansion
- Laguerre polynomial superposition
- Adjustable edge sharpness

**Use Cases:** Uniform illumination, laser processing, material ablation

---

### 8. FrozenWaveBeam
**Beams with arbitrary longitudinal intensity profiles**

```matlab
beam = FrozenWaveBeam(wavelength, w0_aperture, L_pattern, Q_factor, N_modes, target_func)
```

**Parameters:**
- `wavelength`: Wavelength in meters
- `w0_aperture`: Physical aperture radius
- `L_pattern`: Length of pattern segment
- `Q_factor`: Offset parameter (controls transverse spot size)
- `N_modes`: Number of Bessel beams to superpose
- `target_func`: Function handle `@(z)` for desired intensity profile

**Key Features:**
- Fourier-Bessel expansion
- Customizable longitudinal patterns
- Gaussian apodization for physical realizability

**Use Cases:** Custom beam shaping, structured illumination

---

### 9. HelicoConicalBeam
**Helico-conical beams with spiral phase-radial coupling**

```matlab
beam = HelicoConicalBeam(l, K, r0, wavelength, w0_gauss)
```

**Parameters:**
- `l`: Topological charge (OAM)
- `K`: Conical parameter (controls spiral tail winding)
- `r0`: Radial normalization scale
- `wavelength`: Wavelength in meters
- `w0_gauss`: Gaussian apodization radius

**Key Features:**
- Phase: `l*phi + K*(r/r0 - 1)`
- FFT propagation support
- Spiral wavefront structure

**Use Cases:** Specialized beam patterns, optical manipulation

---

### 10. InceGaussianBeam
**Ince-Gaussian beams in elliptic symmetry**

```matlab
beam = InceGaussianBeam(p, m, epsilon, parity, wavelength, w0)
```

**Parameters:**
- `p`: Mode order (≥0)
- `m`: Mode degree (0 ≤ m ≤ p, p-m must be even)
- `epsilon`: Ellipticity parameter (0 = LG limit, ∞ = HG limit)
- `parity`: 'even' or 'odd'
- `wavelength`: Wavelength in meters
- `w0`: Beam waist radius

**Key Features:**
- Numerical eigenvalue solver
- Elliptic coordinate transformation
- Intermediate between LG and HG

**Use Cases:** Elliptical resonator modes, intermediate symmetry beams

---

### 11. MathieuBeam
**Non-diffracting Mathieu beams in elliptic coordinates**

```matlab
beam = MathieuBeam(order, q_param, parity, wavelength, w0_aperture)
```

**Parameters:**
- `order`: Mode order m
- `q_param`: Ellipticity parameter q (≥0)
- `parity`: 'even' or 'odd'
- `wavelength`: Wavelength in meters
- `w0_aperture`: Gaussian aperture scale

**Key Features:**
- Angular and radial Mathieu function evaluation
- Eigenvalue problem solver
- Bessel limit (q=0) handling

**Use Cases:** Non-diffracting elliptic beams, structured light research

---

### 12. PearceyBeam
**Pearcey beams (cusp catastrophe) with auto-focusing**

```matlab
beam = PearceyBeam(wavelength, w0, z_focus)
```

**Parameters:**
- `wavelength`: Wavelength in meters
- `w0`: Transverse scaling parameter
- `z_focus`: Distance where auto-focusing occurs

**Key Features:**
- Numerical Pearcey integral evaluation
- Auto-focusing behavior
- FFT propagation support (static method)

**Use Cases:** Self-healing beams, auto-focusing applications

---

### 13. PerfectVortexBeam
**Perfect vortex with fixed ring radius (OAM-independent)**

```matlab
beam = PerfectVortexBeam(l, r_ring, w0, wavelength)
```

**Parameters:**
- `l`: Topological charge (OAM order)
- `r_ring`: Fixed radius of the bright ring (meters)
- `w0`: Width of the ring (Gaussian thickness)
- `wavelength`: Wavelength in meters

**Key Features:**
- Modified Bessel function implementation
- Ring radius independent of OAM
- Numerical stability with scaled Bessel functions

**Use Cases:** Optical communications with fixed ring size, OAM multiplexing

---

### 14. VectorBeam
**Cylindrical vector beams with spatially varying polarization**

```matlab
beam = VectorBeam(type, wavelength, w0, order)
```

**Parameters:**
- `type`: 'radial', 'azimuthal', or 'hybrid'
- `wavelength`: Wavelength in meters
- `w0`: Beam waist
- `order`: Polarization order P (default: 1)

**Key Features:**
- Spatially varying Jones vectors
- Radial, azimuthal, and hybrid polarization
- Transverse and longitudinal field components

**Use Cases:** Polarization-structured beams, optical manipulation, quantum optics

---

## Analysis Scripts

Each beam type has a corresponding analysis script in the `/analysis/` directory:

| Script | Beam Type | Key Visualizations |
|--------|-----------|-------------------|
| `run_gaussian_analysis.m` | Gaussian | Near/far field, beam caustic |
| `run_lg_analysis.m` | Laguerre-Gaussian | Radial profile, OAM phase, propagation |
| `run_hg_analysis.m` | Hermite-Gaussian | Mode patterns, phase checkerboard |
| `run_bg_analysis.m` | Bessel-Gauss | Diffraction-free zone, longitudinal map |
| `run_airy_analysis.m` | Airy | L-shape intensity, acceleration trajectory |
| `run_bottlebeam_analysis.m` | Bottle | Dark focus visualization |
| `run_flatTop_analysis.m` | Flat-top | Top-hat profile |
| `run_frozen_analysis.m` | Frozen Wave | Custom longitudinal patterns |
| `run_helico_analysis.m` | Helico-conical | Spiral patterns |
| `run_ince_analysis.m` | Ince-Gaussian | Elliptic mode patterns |
| `run_mathieu_analysis.m` | Mathieu | Non-diffracting elliptic patterns |
| `run_pearcey_analysis.m` | Pearcey | Auto-focusing behavior |
| `run_pov_analysis.m` | Perfect Vortex | Fixed-radius ring |
| `run_vector_analysis.m` | Vector | Polarization maps |

**Running Analysis Scripts:**

```matlab
% Navigate to analysis directory
cd analysis

% Run any analysis script
run_gaussian_analysis
run_lg_analysis
% etc.
```

## Common Interface

All beam classes implement a consistent interface:

### Core Methods

```matlab
% Generate complex field
field = beam.generate_beam_field(r, phi, z, ...)
% or for Cartesian beams:
field = beam.generate_beam_field(x, y, z, ...)

% Calculate intensity
intensity = beam.calculate_intensity(r, phi, z, ...)

% Get beam parameters at distance z
params = beam.get_beam_parameters(z)
```

### Common Parameters

Most `generate_beam_field()` methods support:

- `'P_tx_watts'`: Transmit power (default: 1.0 W)
- `'tx_aperture_radius'`: Aperture radius for clipping
- `'beam_tilt_x_rad'`: Beam tilt in x-direction (radians)
- `'beam_tilt_y_rad'`: Beam tilt in y-direction (radians)

### Advanced Parameters (LG, BG, Gaussian)

- `'laser_linewidth_kHz'`: Laser linewidth for phase noise
- `'timing_jitter_ps'`: Timing jitter (picoseconds)
- `'phase_noise_samples'`: Pre-computed phase noise sequence
- `'symbol_time_s'`: Symbol time for noise calculations

## Advanced Features

### Phase Noise Modeling

```matlab
% Generate phase noise sequence
num_symbols = 1000;
symbol_time = 1e-9;  % 1 ns
linewidth = 10.0;     % kHz
phase_noise = beam.generate_phase_noise_sequence(num_symbols, symbol_time, linewidth);

% Use in field generation
field = beam.generate_beam_field(r, phi, z, ...
    'phase_noise_samples', phase_noise(1));
```

### Aperture Loss Calculation

```matlab
% Calculate transmission through aperture
aperture_radius = 25e-3;  % 25 mm
transmission = beam.calculate_tx_aperture_loss(aperture_radius);
loss_dB = -10*log10(transmission);
```

### Beam Overlap

```matlab
% Calculate overlap between two LG beams
beam1 = LaguerreGaussianBeam(0, 1, 1550e-9, 25e-3);
beam2 = LaguerreGaussianBeam(0, 2, 1550e-9, 25e-3);
overlap = beam1.overlap_with(beam2);
```

### Propagation Summary

```matlab
% Get detailed propagation parameters
z_distances = [0, 100, 500, 1000, 2000];
beam.propagation_summary(z_distances);
```

## Examples

### Example 1: Compare Gaussian and LG Beams

```matlab
wavelength = 1550e-9;
w0 = 25e-3;

% Create beams
gaussian = GaussianBeam(wavelength, w0);
lg_beam = LaguerreGaussianBeam(0, 2, wavelength, w0);

% Generate fields
r = linspace(0, 50e-3, 200);
phi = zeros(size(r));
z = 0;

int_gauss = gaussian.calculate_intensity(r, phi, z);
int_lg = lg_beam.calculate_intensity(r, phi, z);

% Plot comparison
figure;
plot(r*1e3, int_gauss, 'b-', 'LineWidth', 2); hold on;
plot(r*1e3, int_lg, 'r-', 'LineWidth', 2);
xlabel('Radius [mm]');
ylabel('Intensity [a.u.]');
legend('Gaussian', 'LG_{0,2}');
grid on;
```

### Example 2: Airy Beam Acceleration

```matlab
wavelength = 532e-9;
x0 = 50e-6;
a = 0.05;
beam = AiryBeam(wavelength, x0, a);

% Calculate trajectory
z = linspace(0, 0.1, 200);
x_traj = beam.trajectory(z);

% Plot parabolic path
plot(z*100, x_traj*1e6, 'r-', 'LineWidth', 2);
xlabel('Propagation Distance z [cm]');
ylabel('Deflection x [μm]');
title('Airy Beam Parabolic Trajectory');
```

### Example 3: Vector Beam Polarization

```matlab
wavelength = 1064e-9;
w0 = 10e-3;
beam = VectorBeam('radial', wavelength, w0, 1);

% Generate vector fields
r_max = 20e-3;
grid_size = 200;
x = linspace(-r_max, r_max, grid_size);
y = linspace(-r_max, r_max, grid_size);
[X, Y] = meshgrid(x, y);
[PHI, R] = cart2pol(X, Y);

[Ex, Ey, ~] = beam.generate_vector_fields(R, PHI, 0);

% Plot intensity
intensity = abs(Ex).^2 + abs(Ey).^2;
imagesc(x*1e3, y*1e3, intensity);
axis xy equal tight;
colormap('hot');
colorbar;
```

## References

### Key Papers

1. **Airy Beams**: Siviloglou & Christodoulides, Opt. Lett. 32, 979 (2007)
2. **Bessel-Gauss**: Gori, Guattari, & Padovani, Opt. Commun. 64, 491 (1987)
3. **Frozen Waves**: Vieira et al., Opt. Lett. 37, 4970 (2012)
4. **Perfect Vortex**: Vaity & Rusch, Opt. Lett. 40, 597 (2015)
5. **Flattened Gaussian**: Gori, Opt. Commun. 107, 335 (1994)
6. **Ince-Gaussian**: Bandres & Gutierrez-Vega, Opt. Lett. 29, 144 (2004)
7. **Pearcey Beams**: Ring et al., Opt. Express 20, 25697 (2012)
8. **Bottle Beams**: Arlt & Padgett, Opt. Lett. 25, 191 (2000)
9. **Helico-conical**: Alonzo et al., Opt. Lett. 30, 3050 (2005)

### General References

- Siegman, A. E. "Lasers" (1986) - Fundamental beam propagation
- Allen et al., "Orbital angular momentum of light" - OAM beams
- Zhan, Q. "Cylindrical vector beams" - Vector beams

## Contributing

Contributions are welcome! Areas for improvement:
- Additional beam types
- Performance optimizations
- Extended documentation
- Unit tests
- Example applications

## License

MIT License - see LICENSE file for details.

## Contact

Contact me at - connect.davuluri@gmail.com

---



