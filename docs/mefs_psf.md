# Scene Simulation (`mefs_psf`)

The `mefs_psf` program is the core simulation engine of the MEFS toolkit. It performs Fourier
optical propagation of point sources through an unobstructed circular pupil, an optional phase-only
coronagraph, and an Integral Field Spectrograph (IFS) consisting of a lenslet array, pinhole mask,
and dispersion optics.

---

## Command Line Options

Run `mefs_psf --help` to display the unified help page.

```
NAME
  mefs_psf - Fourier optical simulation of coronagraphic PSFs and IFS instruments.

USAGE
  mefs_psf [OPTIONS]

OPTIONS
  -h, --help                  Show this help message.
  -s, --single                Run in single on-axis point source mode.
  -S, --scene [<file>]        Run in scene simulation mode. If file is omitted,
                              scans and simulates all 'scene.*.txt' files.
  -c, --coro_order <val>      Coronagraph order (0=unobscured, 2, 4, 6, 8, default: 0).
  -d, --pupil_diam <val>      Pupil diameter in pixels (default: 25.6).
  -g, --grid_size <val>       Pupil grid size (default: 512).
  -l, --lenslet_size <val>    Lenslet size in lambda/D (default: 0.5).
  -n, --lenslet_count <val>   Lenslet grid dimensions (count x count, default: 32).
  -p, --pinhole_diam <val>    Pinhole diameter in lambda/D (0=no pinhole, default: 0.0).
  -f, --ifs_fft_size <val>    IFS FFT grid size (default: 128).
  -m, --ifs_mask <file>       Apply FITS transmission mask to IFS detector output.
```

---

## Simulation Modes

The program runs in two main modes:

### 1. Single-Source Mode (`-s` / `--single`)
Simulates a single, on-axis point source. This mode is useful for generating reference coronagraphic
PSFs, verifying system throughput, or debugging lenslet alignment.

### 2. Scene Simulation Mode (`-S` / `--scene`)
Accepts an astronomical scene file containing a list of point sources. Each line in the scene file
represents a point source with three space-separated values:
```
<x_coordinate_ld>  <y_coordinate_ld>  <flux_intensity>
```

#### Individual Scene Processing
If a file is provided (e.g. `./mefs_psf --scene scene.planet.txt`), `mefs_psf` simulates all sources
in the file, aggregates the output intensities in memory, and writes:
- `psf_im.fits`: Direct focal plane PSF intensity map.
- `ifsraw_im.fits`: Arranged IFS detector plane raw intensity.
- `ifs_im.fits`: Integrated and normalized lenslet flux map.

#### Batch Scene Mode
If `--scene` is run **without arguments** (or simply `./mefs_psf -S`), the program enters batch
mode:
1. **Glob Scanning**: Scans the current working directory for files matching `scene.*.txt` (e.g.
   `scene.star.txt`, `scene.zodi.txt`, `scene.exozodi.txt`, `scene.planet.txt`).
2. **Skip-on-Exist Check**: Checks if the log summary file `scene.<name>.log` exists. If present,
   it skips simulating this scene to save computation time.
3. **Incremental FITS Loading**: When skipping a scene, the program reads the scene's written FITS
   files from disk and adds them to the master sum accumulators, guaranteeing complete and accurate
   combined outputs.
4. **Scenesum Output**: Writes the combined totals of all processed and loaded scenes to:
   - `psf_im.scenesum.fits`
   - `ifsraw_im.scenesum.fits`
   - `ifs_im.scenesum.fits`

---

## MEFS Projection Mode (`-M` / `--mefs`)

The Matched Electric Field Spectrograph (MEFS) mode projects the focal plane complex amplitude
\(b\) (electric field) of every simulated point source onto a reference planet mode \(a\).

The coupled power \(F\) is calculated via the overlap integral:
\[F = \frac{\left|\sum a_i^* b_i\right|^2}{\sum |a_i|^2}\]

When MEFS mode is active:
- **Reference Mode**: The reference planet complex amplitude \(a\) is computed on-axis or
  off-axis by simulating a source of unit flux at coordinates `mefs_x` and `mefs_y`.
- **ASCII Coupling Files**: For each scene, it calculates the coupling efficiency and
  coupled flux for each point, writing them to `mefs.<name>.txt`.
- **Log Summary Report**: Once all scenes are computed or loaded, it aggregates the total planet
  light \(P\) and total background light \(B\) (starlight, zodi, exozodi) to calculate the combined
  spectrograph SNR:
  \[SNR = \frac{P}{\sqrt{P + B}}\]
  The report is printed to the console and written to `mefs.log`.
- **No FITS Files**: MEFS mode bypasses all image creation and does not output any FITS files.

---

## Under-the-Hood Optimizations

- **Parallel Point Loop**: Distributes point sources across CPU cores using OpenMP loops.
  FFTW plans and temporary lenslet pad arrays are allocated as thread-private variables.
- **Resource Reuse**: To maximize throughput, FFT plans, FITS masks, and scaling arrays are
  loaded and pre-allocated once at startup, rather than inside the simulation loops.
- **Header Coordinate Embedding**: The output normalized lenslet flux map (`ifs_im.fits`)
  is written with WCS cards (`CDELT1`, `CDELT2`, `CRPIX1`, `CRPIX2`, `CRVAL1`, `CRVAL2`)
  to specify the physical coordinate system automatically.
