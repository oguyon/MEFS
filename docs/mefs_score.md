# Performance Scorer (`mefs_score`)

The `mefs_score` utility reads the planet flux map and the combined scene sum flux map FITS files,
extracts WCS coordinates to map pixels to physical angular separation in \(\lambda/D\), sorts the
pixels by decreasing planet intensity, and writes a structured table containing cumulative flux
integrals.

This tool is used to evaluate contrast curves and planet detection metrics in the presence of
zodiacal dust, exozodi light, and coronagraphic starlight leaks.

---

## Command Line Options

Run `mefs_score --help` to view options.

```
NAME
  mefs_score - Grade planet detection contrast metric from IFS outputs.

USAGE
  mefs_score [OPTIONS]

OPTIONS
  -p, --planet <file>        Planet flux map FITS file (default: ifs_im.planet.fits).
  -s, --scenesum <file>      Combined scenesum FITS file (default: ifs_im.scenesum.fits).
  -o, --output <file>        Output ASCII file path (default: score.txt).
```

---

## Automatic Coordinate Extraction (WCS)

Unlike other tools that require manual coordinate scales to be passed as inputs, `mefs_score` is
header-driven and reads coordinates directly from the `planet` FITS file:
- **Pixel Scale (`CDELT1`)**: Read from the WCS headers as the angular lenslet size in \(\lambda/D\)
  (defaults to `0.5` \(\lambda/D\) if missing).
- **Reference Center (`CRPIX1`, `CRPIX2`)**: Read from WCS headers as the coordinate origin pixel
  index (defaults to the center of the array if missing).

The physical coordinate of each lenslet pixel index \((lx, ly)\) is calculated using standard WCS
coordinate translation:
\[x_{\lambda/D} = (lx + 1.0 - CRPIX_1) \times CDELT_1\]
\[y_{\lambda/D} = (ly + 1.0 - CRPIX_2) \times CDELT_1\]

---

## Output File Format (`score.txt`)

The generated ASCII file contains a detailed `#` header followed by space-aligned columns for each
pixel in the grid:

| Column | Name | Unit | Description |
| :--- | :--- | :--- | :--- |
| **Col 1** | `Rank index` | integer | 1-based pixel rank sorted by decreasing planet intensity. |
| **Col 2** | `X coordinate` | \(\lambda/D\) | Angular separation along Axis 1. |
| **Col 3** | `Y coordinate` | \(\lambda/D\) | Angular separation along Axis 2. |
| **Col 4** | `Planet intensity` | arbitrary | Flux contribution from the planet alone. |
| **Col 5** | `Total intensity` | arbitrary | Flux from the combined scene sum. |
| **Col 6** | `Difference` | arbitrary | Background light leak (`Total - Planet`). |
| **Col 7** | `Cumulative Planet` | arbitrary | Sum of planet flux up to this point. |
| **Col 8** | `Cumulative Total` | arbitrary | Sum of total flux up to this point. |
| **Col 9** | `Cumulative Diff` | arbitrary | Sum of background leak flux up to this point. |
| **Col 10**| `SNR` | arbitrary | Pixel Signal-to-Noise Ratio: \(P / \sqrt{P+T}\). |
| **Col 11**| `Cumulative SNR` | arbitrary | Combined SNR: \(\sqrt{\sum SNR_i^2}\). |
| **Col 12**| `Row index` | integer | 0-indexed row coordinate in the FITS file. |
| **Col 13**| `Col index` | integer | 0-indexed column coordinate in the FITS file. |
