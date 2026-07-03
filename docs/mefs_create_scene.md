# Scene Generator (`mefs_create_scene`)

The `mefs_create_scene` utility program generates simulated astronomical scenes as lists of
normalized point sources in an ASCII file format. These files are directly consumed by the `mefs_psf`
simulation program.

---

## Command Line Options

Run `mefs_create_scene --help` to view options.

```
NAME
  mefs_create_scene - Generate astronomical scene point source lists.

USAGE
  mefs_create_scene [OPTIONS]

OPTIONS
  -h, --help                  Show this help message.
  -o, --output <file>         Output scene file path (default: scene.txt).
  -t, --type <type>           Scene type: planet, disk, or exozodi.
  -f, --flux <val>            Total integrated flux (default: 1.0).
  -r, --radius <val>          Outer radius for uniform disk in lambda/D (default: 1.0).
  -s, --sampling <val>        Point spacing in lambda/D (default: 0.1).
  -x, --xcoord <val>          Planet X coordinate in lambda/D (default: 1.5).
  -y, --ycoord <val>          Planet Y coordinate in lambda/D (default: 0.0).
  -i, --inclination <val>     Disk inclination angle in degrees (default: 60.0 for exozodi, 0.0 for disk).
  -a, --append                Append generated points to output file instead of overwriting.
  --inner <val>               Exozodi disk inner radius in lambda/D (default: 0.15).
  --outer <val>               Exozodi disk outer radius in lambda/D (default: 3.0).
  --power <val>               Exozodi power law radial index (default: -2.0).
```

---

## Default Benchmark Scenes

When run **without options** (simply `./mefs_create_scene`), the program automatically creates
four standard benchmark files in the current working directory:

1. **`scene.star.txt`**:
   Represents a resolved stellar surface disk.
   - Geometry: Uniform disk, radius = \(0.0075 \lambda/D\).
   - Flux: \(1.0 \times 10^{10}\) integrated.
   - Point Sampling: \(0.001 \lambda/D\) (172 point sources).
2. **`scene.zodi.txt`**:
   Represents local zodiacal dust background.
   - Geometry: Uniform circular disk, radius = \(8.0 \lambda/D\).
   - Flux: \(200.0\) integrated.
   - Point Sampling: \(0.1 \lambda/D\) (20,079 point sources).
3. **`scene.exozodi.txt`**:
   Represents an inclined circumstellar exozodiacal dust disk.
   - Geometry: Elliptical inclined power-law disk, inner radius = \(0.15 \lambda/D\), outer
     radius = \(3.0 \lambda/D\), radial power law index = \(-2.0\).
   - Flux: \(200.0\) integrated.
   - Inclination: \(60.0^{\circ}\) (statistical median inclination).
   - Point Sampling: \(0.1 \lambda/D\) (1,398 point sources).
4. **`scene.planet.txt`**:
   Represents a point source exoplanet.
   - Geometry: Single point source, located at \((1.5, 0.0) \lambda/D\).
   - Flux: \(1.0\).

---

## Exozodiacal Disk Geometry & Inclination

Exozodiacal disks are modeled with a radial power-law surface brightness distribution:
\[I(r) \propto r^{\alpha}\]
where \(\alpha\) is the power-law index (default: \(-2.0\)).

### Inclined Projection Math
To represent an inclined disk in the sky plane:
1. Coordinates are generated on a square grid with step `sampling`.
2. For each point \((x, y)\), we project it into the inclined disk plane coordinate system \((x', y')\):
   - \(x' = x\)
   - \(y' = \frac{y}{\cos i}\)
   where \(i\) is the inclination angle.
3. The physical radius in disk coordinates is \(r = \sqrt{x'^2 + y'^2}\).
4. If \(r\) lies between the inner radius \(r_{in}\) and outer radius \(r_{out}\), the point is
   retained, and its raw intensity is set to:
   \[I(r) = r^{\alpha}\]
5. Finally, the intensities of all retained points are summed and scaled so that the integrated
   total matches the target `--flux` parameter.
