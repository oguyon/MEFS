# Matched Electric Field Spectrograph (MEFS) Toolkit

The Matched Electric Field Spectrograph (MEFS) simulation and analysis toolkit is an end-to-end
numerical package for modeling high-contrast, high-resolution exoplanet spectroscopy instruments
designed for future space telescopes (e.g., Habitable Worlds Observatory - HWO) and ELTs.

Live documentation is available at: **https://oguyon.github.io/MEFS/**

---

## 1. Instrument Concept & Context

MEFS uses a **two-stage starlight suppression strategy** to enable high-resolution spectroscopy
and high-precision wavefront self-calibration:

1. **Stage 1 (Coronagraph)**: Suppresses the bulk of the host starlight ($\sim 10^7$ contrast),
   relaxing requirements on the second stage and narrowing the field of view (FOV) to $\sim 2 \lambda/D$.
2. **Stage 2 (Photonic Integrated Circuit - PIC)**: Utilizes a **Photonic Lantern (PL)** and
   **Photonic Nulling Chip (PNC)** to null residual starlight, decompose electric field modes,
   and route planet light directly to the science spectrograph while sending residual starlight to
   wavefront sensors (WFS) for closed-loop control and self-calibration.

This spatial-spectral modal decomposition distinguishes exoplanet light from starlight speckles
at the science detector, achieving extreme precision ($\sim 50\ \mu\text{as}$ per spectral bin).

---

## 2. Toolkit Structure & Executables

The toolkit compiles into three core executables:

* **[mefs_psf](docs/mefs_psf.md)**:
  Simulates the optical propagation, coronagraph mask, photonic lantern injection, modal nulling,
  and spectrograph coupling.
* **[mefs_create_scene](docs/mefs_create_scene.md)**:
  Generates synthetic astronomical scenes (point sources for stars/planets, disk dust distributions
  for zodiacal and exozodiacal light).
* **[mefs_score](docs/mefs_score.md)**:
  Evaluates spectroscopic signal-to-noise ratios (SNR), coupled/uncoupled flux metrics, and overall
  instrument performance based on simulation logs.

---

## 3. Quick Start & Building

### Prerequisites
- CMake ($\ge 3.10$)
- GCC or Clang (supporting C11 and OpenMP)
- CFITSIO
- FFTW3

### Build Instructions
```bash
mkdir -p build
cd build
cmake ..
make -j$(nproc)
```

### Running Simulations
1. Generate the standard benchmark scene files:
   ```bash
   ./build/mefs_create_scene --benchmark
   ```
2. Run the PSF simulation using the planet scene:
   ```bash
   ./build/mefs_psf --scene scene.planet.txt --mefs_x 1.0
   ```
3. Evaluate the output SNR:
   ```bash
   ./build/mefs_score --scene scene.planet.txt
   ```

---

## 4. References & Project Status

MEFS technology is actively prototyped and demonstrated on-sky via the **SCExAO** instrument at the
Subaru Telescope:
- **GLINT (Photonic Nulling)**: On-sky demonstration (May 2025) showing active phase control, WFS/C
  loop-closing, and $>80\%$ through-chip transmission.
- **FIRST / Photonic Lantern**: Demonstrated spatial-spectral modal decomposition and sub-diffraction
  limited astronomy (*Yoo Jung Kim et al. 2025 ApJL*).
