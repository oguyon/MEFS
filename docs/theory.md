# Physical & Optical Theory

This page outlines the physical, optical, and statistical concepts implemented in the MEFS toolkit,
incorporating the two-stage starlight suppression strategy, photonic nulling, and wavefront
self-calibration.

---

## 1. Lyot Coronagraphs & Stage 1 Suppression

A coronagraph is an optical instrument designed to suppress the coherent diffraction of light from a
bright central starlight to reveal faint companion planets or circumstellar dust disks.

In a **Phase Mask Lyot Coronagraph**, a phase shift is applied in the focal plane to diffract
starlight outside of the pupil boundary in a downstream pupil plane, where a Lyot stop blocks it.

In the MEFS instrument architecture:
* **Stage 1 (Coronagraph)**: Suppresses the bulk of the host starlight ($\sim 10^7$ contrast
  reduction).
* This relaxes the Stage 2 contrast requirement to only $\sim 10^2 - 10^3$.
* It also reduces the required Field of View (FOV) to $\sim 2\ \lambda/D$, reducing the number of
  photonic channels by $\sim 25\times$.

---

## 2. Stage 2: Photonic Integrated Circuit (PIC)

Stage 2 performs fine nulling and modal decomposition using a **Photonic Lantern (PL)** and
**Photonic Nulling Chip (PNC)**:

1. **Photonic Lantern (PL)**: A multi-mode interface is injected with the focal plane complex
   amplitude. It decomposes the spatial electric field into a discrete set of single-mode fiber
   (SMF) cores.
2. **Photonic Nulling Chip (PNC)**: Applies precise phase shifts (using couplers and phase
   shifters) to null starlight.
3. **Science Routing**: Exoplanet light is routed directly to the science spectrograph, while
   residual starlight leakage is sent to wavefront sensors (WFS) for closed-loop tracking and
   correction.

---

## 3. Spatial-Spectral Modal Decomposition & Self-Calibration

Instead of trying to achieve perfect starlight suppression purely through optical alignment, MEFS
measures the starlight leakage directly at the science detector.

The coupled power \(F\) is calculated via the overlap integral:
\[F = \frac{\left|\sum a_i^* b_i\right|^2}{\sum |a_i|^2}\]
where \(b\) is the complex electric field amplitude of the scene, and \(a\) is the reference mode.

### Wavefront Sensing & Precision
By decomposing the electric field into spatial-spectral modes, WFS/C tracking loops can distinguish
between actual exoplanet signals and starlight speckles caused by wavefront errors.
* This enables high-precision self-calibration, achieving up to **$50\ \mu\text{as}$ precision per
  spectral bin** on-sky.
* It relaxes sub-nanometer wavefront control stability requirements compared to classical
  coronagraphy.

---

## 4. Fractional Coronagraph Mode Subtraction

The toolkit models starlight leakage modes and allows subtracting them with fractional weights,
interpolating continuously between 0 and 6 modes using `coronagraph_param` ($P_c \in [0.0, 3.0]$):

### Mode Weights
Let $w_m$ represent the subtraction weight of the $m$-th leakage mode (for $m = 0, \dots, 5$):
* **$P_c \le 0$**: No modes are subtracted ($w_m = 0$ for all $m$).
* **$0 < P_c \le 1$**: Subtracted mode 0 (starlight mode) with weight $w_0 = P_c$.
* **$1 < P_c \le 2$**: Mode 0 is fully subtracted ($w_0 = 1.0$); modes 1 and 2 are subtracted with
  fractional weight $w_1 = w_2 = P_c - 1.0$.
* **$2 < P_c \le 3$**: Modes 0, 1, 2 are fully subtracted ($w_0 = w_1 = w_2 = 1.0$); modes 3, 4, 5
  are subtracted with fractional weight $w_3 = w_4 = w_5 = \min(P_c - 2.0, 1.0)$.

This continuous formulation enables detailed starlight suppression optimization studies.

---

## 5. Circumstellar Dust Disks & Inclination

Circumstellar dust is a major source of background noise in high-contrast imaging. The toolkit
models local zodiacal dust and inclined exozodiacal disks:

### Exozodiacal Dust Model
Modeled as an inclined power-law disk with radial index $\alpha$ (default: $-2.0$).

#### Statistical Median Inclination
If exozodiacal disks are oriented randomly in 3D space, their inclination angles \(i\) follow a
probability distribution:
\[P(i) = \sin i\]
The cumulative distribution function (CDF) representing the probability that a disk has an
inclination angle less than or equal to \(i\) is:
\[F(i) = \int_{0}^{i} \sin \theta \, d\theta = 1 - \cos i\]
The median inclination angle \(i_{med}\) (where half of the randomly oriented disks are less
inclined) is found by setting the CDF to \(0.5\):
\[1 - \cos(i_{med}) = 0.5 \implies \cos(i_{med}) = 0.5 \implies i_{med} = 60.0^{\circ}\]

This statistical median inclination of **\(60.0^{\circ}\)** is adopted as the physical default
for exozodi simulations.

---

## 6. On-Sky Prototyping & SCExAO Status

MEFS technology is actively developed and validated on-sky via the **SCExAO** instrument on the
Subaru Telescope:
* **GLINT (Photonic Nulling)**: Successfully demonstrated active phase control, closed-loop
  WFS/C operations, and $>80\%$ through-chip transmission using phase shifters and tricouplers
  (May 2025).
* **FIRST / Photonic Lantern**: Demonstrated spatial-spectral modal decomposition and
  self-calibration, achieving sub-diffraction-limited measurements (*Yoo Jung Kim et al. 2025 ApJL*).
