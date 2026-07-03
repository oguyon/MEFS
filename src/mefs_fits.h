#ifndef MEFS_FITS_H
#define MEFS_FITS_H

#include <fftw3.h>

/**
 * @brief Read a 2D FITS image into a double array.
 *
 * @param filename File path of the FITS image.
 * @param out_rows Pointer to store the number of rows read.
 * @param out_cols Pointer to store the number of columns read.
 * @return double* Pointer to allocated 2D array of doubles, or NULL on failure.
 */
double *read_fits_2d(
    const char *filename,
    int        *out_rows,
    int        *out_cols);

/**
 * @brief Read a 2D or 3D FITS image as a complex amplitude array.
 *
 * The input is expected to contain Plane 0: Amplitude and Plane 1: Phase (optional, in radians).
 *
 * @param filename File path of the FITS image.
 * @param out_rows Pointer to store the number of rows read.
 * @param out_cols Pointer to store the number of columns read.
 * @return fftw_complex* Pointer to allocated complex array, or NULL on failure.
 */
fftw_complex *read_fits_complex(
    const char *filename,
    int        *out_rows,
    int        *out_cols);

/**
 * @brief Write a 2D double array as a 2D FITS image.
 *
 * @param filename File path of the output FITS image.
 * @param data Flat 2D array of double data.
 * @param rows Number of rows in the grid.
 * @param cols Number of columns in the grid.
 */
void write_fits_2d(
    const char   *filename,
    const double *data,
    int           rows,
    int           cols);

/**
 * @brief Write a 2D double array with WCS coordinates as a FITS image.
 *
 * @param filename File path of the output FITS image.
 * @param data Flat 2D array of double data.
 * @param rows Number of rows in the grid.
 * @param cols Number of columns in the grid.
 * @param cdelt Pixel scale delta value for WCS header.
 * @param crpix Reference pixel coordinate center for WCS header.
 */
void write_fits_2d_wcs(
    const char   *filename,
    const double *data,
    int           rows,
    int           cols,
    double        cdelt,
    double        crpix);

/**
 * @brief Write a complex array as a 3D FITS image (Plane 0: Amplitude, Plane 1: Phase).
 *
 * @param filename File path of the output FITS image.
 * @param data Flat complex array.
 * @param rows Number of rows in the grid.
 * @param cols Number of columns in the grid.
 */
void write_fits_complex(
    const char         *filename,
    const fftw_complex *data,
    int                 rows,
    int                 cols);

#endif /* MEFS_FITS_H */
