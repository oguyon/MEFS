# Physical & Optical Theory

This page outlines the physical, optical, and statistical concepts implemented in the MEFS toolkit.

---

## Lyot Coronagraphs

A coronagraph is an optical instrument designed to suppress the coherent diffraction of light from a
bright central star to reveal faint companion planets or circumstellar dust disks.

In a **Phase Mask Lyot Coronagraph**, a phase shift is applied in the focal plane to diffract starlight
outside of the pupil boundary in a downstream pupil plane, where a Lyot stop blocks the starlight.

The toolkit models phase masks of orders 2, 4, 6, and 8. The transmission in the focal plane is represented
conceptually as:
\[T(\theta) = 1 - \exp\left( - \left(\frac{\theta}{\theta_0}\right)^{2k} \right)\]
where \(2k\) represents the coronagraph order (e.g. 2, 4, 6, or 8). Higher order coronagraphs offer better
tolerance to stellar angular size and tip-tilt jitter but result in larger inner working angles (IWA).

---

## Integral Field Spectrograph (IFS) Simulation

An IFS slices a 2D spatial field and disperses each spatial element (lenslet) into a spectrum on a 2D
detector. In the monochromatic case modeled by this toolkit, the simulation is as follows:

1. **Lenslet Slicing**: The focal plane complex amplitude field is tiled into a grid of \(M \times M\)
   square lenslets of physical size `lenslet_size_ld` \(\lambda/D\).
2. **Pinhole Transmission**: Inside each lenslet, the field is optionally spatial-filtered by a
   circular pinhole mask (representing a micro-lens focal plane array mask).
3. **Propagation to Detector**: The light from each lenslet is propagated to the detector plane via
   a Fast Fourier Transform (FFT) of size \(P \times P\) (`ifs_fft_size`). The output intensity is
   written as the raw IFS detector image.

---

## Circumstellar Dust Disks & Inclination

Circumstellar dust is a major source of background noise in high-contrast imaging. The toolkit models
two types of dust:

### 1. Local Zodiacal Dust
Modeled as a uniform background disk.

### 2. Exozodiacal Dust
Modeled as an inclined power-law disk.

#### Statistical Median Inclination
If exozodiacal disks are oriented randomly in 3D space, their inclination angles \(i\) follow a
probability distribution:
\[P(i) = \sin i\]
The cumulative distribution function (CDF) representing the probability that a disk has an inclination
angle less than or equal to \(i\) is:
\[F(i) = \int_{0}^{i} \sin \theta \, d\theta = 1 - \cos i\]
The median inclination angle \(i_{med}\) (where half of the randomly oriented disks are less inclined and
half are more inclined) is found by setting the CDF to \(0.5\):
\[1 - \cos(i_{med}) = 0.5 \implies \cos(i_{med}) = 0.5 \implies i_{med} = 60.0^{\circ}\]

This statistical median inclination of **\(60.0^{\circ}\)** is adopted as the physical default for exozodi
simulations in `mefs_create_scene`.
