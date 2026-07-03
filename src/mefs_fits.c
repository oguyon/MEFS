#include "mefs_fits.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fitsio.h>
#include <fftw3.h>

/**
 * @brief Read a 2D FITS image into a double array.
 */
double *read_fits_2d(
    const char *filename,
    int        *out_rows,
    int        *out_cols)
{
    fitsfile *fptr = NULL;
    int status = 0;

    fits_open_file(&fptr, filename, READONLY, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error opening FITS file %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
        return NULL;
    }

    int naxis = 0;
    fits_get_img_dim(fptr, &naxis, &status);

    long naxes[2] = {1, 1};
    fits_get_img_size(fptr, 2, naxes, &status);

    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error reading image dimensions for %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return NULL;
    }

    *out_rows = naxes[1];
    *out_cols = naxes[0];

    double *data = malloc(naxes[0] * naxes[1] * sizeof(double));
    if (data == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Failed to allocate memory for image data.\n",
                __func__, __LINE__);
        fits_close_file(fptr, &status);
        return NULL;
    }

    long fpixel[2] = {1, 1};
    fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0] * naxes[1], NULL, data, NULL, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error reading pixels from FITS: ",
                __func__, __LINE__);
        fits_report_error(stderr, status);
        free(data);
        fits_close_file(fptr, &status);
        return NULL;
    }

    fits_close_file(fptr, &status);
    return data;
} // double *read_fits_2d

/**
 * @brief Read a 2D or 3D FITS image as a complex amplitude array.
 */
fftw_complex *read_fits_complex(
    const char *filename,
    int        *out_rows,
    int        *out_cols)
{
    fitsfile *fptr = NULL;
    int status = 0;

    fits_open_file(&fptr, filename, READONLY, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error opening FITS file %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
        return NULL;
    }

    int naxis = 0;
    fits_get_img_dim(fptr, &naxis, &status);

    long naxes[3] = {1, 1, 1};
    fits_get_img_size(fptr, 3, naxes, &status);

    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error reading image dimensions for %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return NULL;
    }

    int cols = naxes[0];
    int rows = naxes[1];
    int planes = (naxis == 3) ? naxes[2] : 1;

    if (naxis < 2 || naxis > 3)
    {
        fprintf(stderr, "ERROR: [%s:%d] FITS image must be 2D or 3D (got %d dimensions).\n",
                __func__, __LINE__, naxis);
        fits_close_file(fptr, &status);
        return NULL;
    }

    if (naxis == 3 && planes != 1 && planes != 2)
    {
        fprintf(stderr, "ERROR: [%s:%d] 3D FITS image must have 1 or 2 planes (got %d).\n",
                __func__, __LINE__, planes);
        fits_close_file(fptr, &status);
        return NULL;
    }

    *out_rows = rows;
    *out_cols = cols;

    fftw_complex *data = fftw_alloc_complex(rows * cols);
    if (data == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Failed to allocate memory for image data.\n",
                __func__, __LINE__);
        fits_close_file(fptr, &status);
        return NULL;
    }

    double *temp_amp = malloc(rows * cols * sizeof(double));
    double *temp_phase = calloc(rows * cols, sizeof(double));

    if (temp_amp == NULL || temp_phase == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Failed to allocate temporary buffers.\n",
                __func__, __LINE__);
        if (temp_amp) free(temp_amp);
        if (temp_phase) free(temp_phase);
        fftw_free(data);
        fits_close_file(fptr, &status);
        return NULL;
    }

    /* Read Amplitude plane (plane 0) */
    long fpixel[3] = {1, 1, 1};
    fits_read_pix(fptr, TDOUBLE, fpixel, rows * cols, NULL, temp_amp, NULL, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error reading Amplitude plane from FITS: ",
                __func__, __LINE__);
        fits_report_error(stderr, status);
        free(temp_amp);
        free(temp_phase);
        fftw_free(data);
        fits_close_file(fptr, &status);
        return NULL;
    }

    /* Read Phase plane (plane 1) if complex */
    if (naxis == 3 && planes == 2)
    {
        fpixel[2] = 2;
        fits_read_pix(fptr, TDOUBLE, fpixel, rows * cols, NULL, temp_phase, NULL, &status);
        if (status)
        {
            fprintf(stderr, "WARNING: [%s:%d] Error reading Phase plane (setting to 0): ",
                    __func__, __LINE__);
            fits_report_error(stderr, status);
            status = 0; /* Clear error status */
        }
    }

    /* Convert Amplitude and Phase to Real and Imaginary for FFTW */
    for (int ii = 0; ii < rows * cols; ii++)
    {
        double amp = temp_amp[ii];
        double phase = temp_phase[ii];
        data[ii][0] = amp * cos(phase);
        data[ii][1] = amp * sin(phase);
    }

    free(temp_amp);
    free(temp_phase);
    fits_close_file(fptr, &status);
    return data;
} // fftw_complex *read_fits_complex

/**
 * @brief Write a 2D double array as a 2D FITS image.
 */
void write_fits_2d(
    const char   *filename,
    const double *data,
    int           rows,
    int           cols)
{
    fitsfile *fptr = NULL;
    int status = 0;

    remove(filename); /* Delete if exists to overwrite */

    fits_create_file(&fptr, filename, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error creating FITS file %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
        return;
    }

    int bitpix = DOUBLE_IMG;
    int naxis = 2;
    long naxes[2] = {cols, rows};
    fits_create_img(fptr, bitpix, naxis, naxes, &status);

    long fpixel[2] = {1, 1};
    fits_write_pix(fptr, TDOUBLE, fpixel, rows * cols, (void *)data, &status);

    fits_close_file(fptr, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error closing FITS file %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
    }
} // void write_fits_2d

/**
 * @brief Write a 2D double array with WCS coordinates as a FITS image.
 */
void write_fits_2d_wcs(
    const char   *filename,
    const double *data,
    int           rows,
    int           cols,
    double        cdelt,
    double        crpix)
{
    fitsfile *fptr = NULL;
    int status = 0;

    remove(filename); /* Delete if exists to overwrite */

    fits_create_file(&fptr, filename, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error creating FITS file %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
        return;
    }

    int bitpix = DOUBLE_IMG;
    int naxis = 2;
    long naxes[2] = {cols, rows};
    fits_create_img(fptr, bitpix, naxis, naxes, &status);

    long fpixel[2] = {1, 1};
    fits_write_pix(fptr, TDOUBLE, fpixel, rows * cols, (void *)data, &status);

    /* Write WCS keywords */
    fits_write_key(fptr, TDOUBLE, "CDELT1", &cdelt, "Pixel scale Axis 1", &status);
    fits_write_key(fptr, TDOUBLE, "CDELT2", &cdelt, "Pixel scale Axis 2", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX1", &crpix, "Reference pixel Axis 1", &status);
    fits_write_key(fptr, TDOUBLE, "CRPIX2", &crpix, "Reference pixel Axis 2", &status);
    double crval = 0.0;
    fits_write_key(fptr, TDOUBLE, "CRVAL1", &crval, "Reference value Axis 1", &status);
    fits_write_key(fptr, TDOUBLE, "CRVAL2", &crval, "Reference value Axis 2", &status);

    fits_close_file(fptr, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error closing FITS file %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
    }
} // void write_fits_2d_wcs

/**
 * @brief Write a complex array as a 3D FITS image.
 */
void write_fits_complex(
    const char         *filename,
    const fftw_complex *data,
    int                 rows,
    int                 cols)
{
    fitsfile *fptr = NULL;
    int status = 0;

    remove(filename); /* Delete if exists to overwrite */

    fits_create_file(&fptr, filename, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error creating complex FITS file %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
        return;
    }

    int bitpix = DOUBLE_IMG;
    int naxis = 3;
    long naxes[3] = {cols, rows, 2};
    fits_create_img(fptr, bitpix, naxis, naxes, &status);

    double *temp = malloc(rows * cols * sizeof(double));
    if (temp == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Failed to allocate memory for FITS serialization.\n",
                __func__, __LINE__);
        fits_close_file(fptr, &status);
        return;
    }

    /* Compute and Write Amplitude plane */
    for (int ii = 0; ii < rows * cols; ii++)
    {
        double r = data[ii][0];
        double im = data[ii][1];
        temp[ii] = sqrt(r * r + im * im);
    }
    long fpixel[3] = {1, 1, 1};
    fits_write_pix(fptr, TDOUBLE, fpixel, rows * cols, temp, &status);

    /* Compute and Write Phase plane (in radians) */
    for (int ii = 0; ii < rows * cols; ii++)
    {
        double r = data[ii][0];
        double im = data[ii][1];
        temp[ii] = atan2(im, r);
    }
    fpixel[2] = 2;
    fits_write_pix(fptr, TDOUBLE, fpixel, rows * cols, temp, &status);

    free(temp);
    fits_close_file(fptr, &status);
    if (status)
    {
        fprintf(stderr, "ERROR: [%s:%d] Error closing complex FITS file %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
    }
} // void write_fits_complex
