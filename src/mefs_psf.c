#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <fitsio.h>
#include <fftw3.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <glob.h>
#include <sys/stat.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Long option character codes for getopt_long */
#define OPT_IFS_COMPLEX   1000
#define OPT_IFS_INTENSITY 1001
#define OPT_IFS_MASK      1002
#define OPT_IFS_FLUX      1003
#define OPT_MEFS_X        1004
#define OPT_MEFS_Y        1005

const char *c_reset   = "";
const char *c_bold    = "";
const char *c_cyan    = "";
const char *c_b_green = "";
const char *c_green   = "";
const char *c_magenta = "";
const char *c_yellow  = "";
const char *c_b_red   = "";
const char *c_grey    = "";

void init_colors(void)
{
    if (getenv("NO_COLOR") == NULL)
    {
        c_reset   = "\x1b[0m";
        c_bold    = "\x1b[1m";
        c_cyan    = "\x1b[1;36m";
        c_b_green = "\x1b[1;32m";
        c_green   = "\x1b[32m";
        c_magenta = "\x1b[35m";
        c_yellow  = "\x1b[33m";
        c_b_red   = "\x1b[1;31m";
        c_grey    = "\x1b[90m";
    }
}

/* Print a progress bar for scene point calculations */
void print_progress_bar(
    int count,
    int total)
{
    int width = 40;
    double progress = (double)(count + 1) / total;
    int filled = (int)(progress * width);

    printf("\r%sProgress:%s %s[%s", c_cyan, c_reset, c_grey, c_reset);
    for (int ii = 0; ii < width; ii++)
    {
        if (ii < filled)
        {
            printf("%s=%s", c_green, c_reset);
        }
        else if (ii == filled)
        {
            printf("%s>%s", c_green, c_reset);
        }
        else
        {
            printf(" ");
        }
    }
    printf("%s]%s %d/%d (%.1f%%)", c_grey, c_reset, count + 1, total, progress * 100.0);
    fflush(stdout);
}

/* Print usage instructions */
void print_usage(
    const char *prog_name)
{
    printf("%sNAME%s\n", c_cyan, c_reset);
    printf("  %s%s%s - Compute Point Spread Functions (PSFs) and "
           "Integral Field Spectrograph (IFS) arrays.\n\n",
           c_b_green, prog_name, c_reset);

    printf("%sUSAGE%s\n", c_cyan, c_reset);
    printf("  %s%s%s [OPTIONS]\n\n", c_b_green, prog_name, c_reset);

    printf("%sDESCRIPTION%s\n", c_cyan, c_reset);
    printf("  Generates telescope pupil plane functions, applies wavefront propagation (FFT)\n");
    printf("  to compute PSFs, and optionally simulates an IFS lenslet array on the PSF\n");
    printf("  focal plane. Supports off-axis sources, multi-order coronagraph suppression,\n");
    printf("  and point source scene list simulations.\n\n");

    printf("%sOPTIONS%s\n", c_cyan, c_reset);
    printf("  %s-i%s %s<file>%s                Input pupil FITS file (plane 0: Amp, "
           "plane 1: Phase).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("                       If omitted, a circular pupil is generated.\n");
    printf("  %s-o%s %s<file>%s                Output complex PSF FITS file "
           "(default: psf_complex.fits).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-p%s %s<file>%s                Output saved pupil wavefront FITS file "
           "(default: pupil.fits).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-t%s %s<file>%s                Output intensity PSF FITS file "
           "(default: psf_intensity.fits).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-n%s %s<size>%s                Grid size for generated circular pupil "
           "(default: 512).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-d%s %s<diam>%s                Diameter in pixels for generated circular "
           "pupil (default: size/20.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-x, --source_x%s %s<val>%s     Off-axis source X position in units of lambda/D "
           "(default: 0.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-y, --source_y%s %s<val>%s     Off-axis source Y position in units of lambda/D "
           "(default: 0.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-c, --coronagraph%s %s<ord>%s  Enable coronagraph subtraction of order 1, "
           "2, or 3 (default: none).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-s, --ifs%s                    Enable the Integral Field Spectrograph (IFS) "
           "simulation.\n",
           c_green, c_reset);
    printf("  %s-l, --lenslet_size%s %s<val>%s  Lenslet size in lambda/D (default: 0.5).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-m, --lenslet_count%s %s<num>%s Lenslet array grid dimensions (default: 32).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-P, --ifs_fft_size%s %s<size>%s  IFS far-field FFT size per lenslet "
           "(default: 32).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-k, --pinhole_diam%s %s<val>%s  Circular pinhole diameter in lambda/dl "
           "(default: 2.0).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s--ifs_mask%s %s<file>%s        Path to 2D FITS file containing transmission mask "
           "overlay.\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s--ifs_complex%s %s<file>%s     Output complex IFS FITS file "
           "(default: ifs_complex.fits).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s--ifs_intensity%s %s<file>%s   Output intensity IFS FITS file "
           "(default: ifs_intensity.fits).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s--ifs_flux%s %s<file>%s        Output 2D lenslet integrated flux FITS file "
           "(default: ifs_flux.fits).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-M, --mefs%s                   Enable MEFS projection mode (matched filter).\n",
           c_green, c_reset);
    printf("  %s--mefs_x%s %s<val>%s            MEFS projection reference mode X coordinate\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("                       in lambda/D (default: 1.0).\n");
    printf("  %s--mefs_y%s %s<val>%s            MEFS projection reference mode Y coordinate\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("                       in lambda/D (default: 0.0).\n");
    printf("  %s-S, --scene%s %s<file>%s        Input ASCII scene definition file "
           "(lines: x y flux; default: none).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-h, --help%s                   Show this help message.\n\n",
           c_green, c_reset);

    printf("%sEXAMPLES%s\n", c_cyan, c_reset);
    printf("  1. Run default IFS simulation:\n");
    printf("     %s$%s %s%s -s%s\n",
           c_grey, c_reset, c_b_green, prog_name, c_reset);
    printf("  2. Compute 3rd order coronagraph PSF off-axis:\n");
    printf("     %s$%s %s%s -c 3 -x 2.0 -y 0.5%s\n",
           c_grey, c_reset, c_b_green, prog_name, c_reset);
    printf("  3. Run multi-source scene simulation:\n");
    printf("     %s$%s %s%s -S scene.txt -s%s\n",
           c_grey, c_reset, c_b_green, prog_name, c_reset);
    printf("  4. Run scene list in MEFS projection mode:\n");
    printf("     %s$%s %s%s -S scene.txt -s -M%s\n\n",
           c_grey, c_reset, c_b_green, prog_name, c_reset);

    printf("%sCOLOR MODE%s\n", c_cyan, c_reset);
    printf("  Colors: %s%s%s. Set %sNO_COLOR%s env var to disable "
           "(see https://no-color.org).\n",
           c_yellow, (getenv("NO_COLOR") == NULL) ? "ENABLED" : "DISABLED", c_reset,
           c_bold, c_reset);
}

/* 2D fftshift/ifftshift operation (out-of-place)
   Moves zero-frequency component to center of spectrum, or vice versa
   (identical for even sizes). */
void fftshift(
    const fftw_complex *in,
    fftw_complex       *out,
    int                 rows,
    int                 cols)
{
    int r_shift = rows / 2;
    int c_shift = cols / 2;

    for (int ii = 0; ii < rows; ii++)
    {
        int target_r = (ii + r_shift) % rows;
        for (int jj = 0; jj < cols; jj++)
        {
            int target_c = (jj + c_shift) % cols;
            int src_idx = ii * cols + jj;
            int dst_idx = target_r * cols + target_c;

            out[dst_idx][0] = in[src_idx][0];
            out[dst_idx][1] = in[src_idx][1];
        }
    }
}

/* Solve linear system G * c = b using Gaussian elimination with partial pivoting.
   N is system size. G is an N x N matrix, b is an N-element vector.
   On return, b contains the solution c. G is modified in place. */
int solve_linear_system(
    double *G,
    double *b,
    int     N)
{
    for (int ii = 0; ii < N; ii++)
    {
        /* Find pivot row */
        int pivot = ii;
        double max_val = fabs(G[ii * N + ii]);
        for (int row = ii + 1; row < N; row++)
        {
            double val = fabs(G[row * N + ii]);
            if (val > max_val)
            {
                max_val = val;
                pivot = row;
            }
        }

        /* Swap rows in G and b */
        if (pivot != ii)
        {
            for (int col = 0; col < N; col++)
            {
                double temp = G[ii * N + col];
                G[ii * N + col] = G[pivot * N + col];
                G[pivot * N + col] = temp;
            }
            double temp = b[ii];
            b[ii] = b[pivot];
            b[pivot] = temp;
        }

        /* Check singularity */
        if (fabs(G[ii * N + ii]) < 1e-12)
        {
            return -1; /* Singular matrix */
        }

        /* Eliminate columns below diagonal */
        for (int row = ii + 1; row < N; row++)
        {
            double factor = G[row * N + ii] / G[ii * N + ii];
            for (int col = ii; col < N; col++)
            {
                G[row * N + col] -= factor * G[ii * N + col];
            }
            b[row] -= factor * b[ii];
        }
    }

    /* Back-substitution */
    for (int ii = N - 1; ii >= 0; ii--)
    {
        for (int col = ii + 1; col < N; col++)
        {
            b[ii] -= G[ii * N + col] * b[col];
        }
        b[ii] /= G[ii * N + ii];
    }

    return 0;
}

void get_coronagraph_weights(
    double coronagraph_param,
    double *weights)
{
    for (int m = 0; m < 6; m++)
    {
        weights[m] = 0.0;
    }

    if (coronagraph_param <= 0.0)
    {
        return;
    }

    if (coronagraph_param <= 1.0)
    {
        weights[0] = coronagraph_param;
    }
    else if (coronagraph_param <= 2.0)
    {
        double x = coronagraph_param - 1.0;
        weights[0] = 1.0;
        weights[1] = x;
        weights[2] = x;
    }
    else
    {
        double x = coronagraph_param - 2.0;
        if (x > 1.0)
        {
            x = 1.0;
        }
        weights[0] = 1.0;
        weights[1] = 1.0;
        weights[2] = 1.0;
        weights[3] = x;
        weights[4] = x;
        weights[5] = x;
    }
}

/* Generate a 2D complex unobstructed circular pupil with sub-pixel anti-aliasing */
fftw_complex *generate_circular_pupil(
    int    rows,
    int    cols,
    double diameter)
{
    fftw_complex *pupil = fftw_alloc_complex(rows * cols);
    if (pupil == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Failed to allocate memory for circular pupil.\n",
                __func__, __LINE__);
        return NULL;
    }

    double cx = (cols - 1) / 2.0;
    double cy = (rows - 1) / 2.0;
    double R = diameter / 2.0;

    for (int ii = 0; ii < rows; ii++)
    {
        double dy = ii - cy;
        for (int jj = 0; jj < cols; jj++)
        {
            double dx = jj - cx;
            double dist = sqrt(dx * dx + dy * dy);
            int idx = ii * cols + jj;

            if (dist < R - 0.7072)
            {
                pupil[idx][0] = 1.0;
            }
            else if (dist > R + 0.7072)
            {
                pupil[idx][0] = 0.0;
            }
            else
            {
                /* Edge pixel: Integrate sub-pixels (20x20 grid) */
                int count = 0;
                int K = 20;

                for (int ky = 0; ky < K; ky++)
                {
                    double sy = ii - 0.5 + (ky + 0.5) / K - cy;
                    for (int kx = 0; kx < K; kx++)
                    {
                        double sx = jj - 0.5 + (kx + 0.5) / K - cx;
                        if (sx * sx + sy * sy <= R * R)
                        {
                            count++;
                        }
                    }
                }
                pupil[idx][0] = (double)count / (K * K);
            }
            pupil[idx][1] = 0.0;
        }
    }
    return pupil;
}

/* Estimate the pupil diameter in pixels from its amplitude profile */
double estimate_pupil_diameter(
    const fftw_complex *pupil,
    int                 rows,
    int                 cols)
{
    int min_r = rows;
    int max_r = -1;
    int min_c = cols;
    int max_c = -1;

    for (int ii = 0; ii < rows; ii++)
    {
        for (int jj = 0; jj < cols; jj++)
        {
            int idx = ii * cols + jj;
            double r = pupil[idx][0];
            double im = pupil[idx][1];
            double amp = sqrt(r * r + im * im);
            if (amp > 0.05) /* 5% threshold */
            {
                if (ii < min_r) min_r = ii;
                if (ii > max_r) max_r = ii;
                if (jj < min_c) min_c = jj;
                if (jj > max_c) max_c = jj;
            }
        }
    }

    if (max_r == -1 || max_c == -1)
    {
        return -1.0;
    }

    double diam_r = max_r - min_r + 1;
    double diam_c = max_c - min_c + 1;
    return (diam_r + diam_c) / 2.0;
}

/* Read 2D double array from FITS image (useful for masks) */
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
}

/* Read 2D/3D FITS into complex fftw_complex array.
   Input format is expected as Plane 0: Amplitude, Plane 1: Phase (in radians). */
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
}

/* Write 2D double array as FITS image */
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
}

/* Write 2D double array with WCS coordinates as FITS image */
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
}

/* Write complex fftw_complex array as a 3D FITS image (Plane 0: Amplitude, Plane 1: Phase) */
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
}

/* Core IFS computation step for a single point source, accumulating into the raw intensity array */
void accumulate_ifs_source(
    const fftw_complex *psf,
    int                 grid_size,
    int                 L,
    int                 M,
    int                 P,
    fftw_complex       *lenslet_padded,
    fftw_complex       *lenslet_padded_shifted,
    fftw_complex       *farfield_shifted,
    fftw_complex       *farfield,
    fftw_plan           ifs_plan,
    double              pinhole_diam_ldl,
    double             *ifsraw_im_accum)
{
    int out_size = M * P;
    int total_lenslet_grid_width = M * L;
    int row_start = (grid_size - total_lenslet_grid_width) / 2;
    int col_start = (grid_size - total_lenslet_grid_width) / 2;

    for (int ly = 0; ly < M; ly++)
    {
        int r0 = row_start + ly * L;
        for (int lx = 0; lx < M; lx++)
        {
            int c0 = col_start + lx * L;
            int pad_start_y = (P - L) / 2;
            int pad_start_x = (P - L) / 2;

            for (int ii = 0; ii < P * P; ii++)
            {
                lenslet_padded[ii][0] = 0.0;
                lenslet_padded[ii][1] = 0.0;
            }

            for (int rr = 0; rr < L; rr++)
            {
                int psf_r = r0 + rr;
                if (psf_r < 0 || psf_r >= grid_size) continue;

                for (int cc = 0; cc < L; cc++)
                {
                    int psf_c = c0 + cc;
                    if (psf_c < 0 || psf_c >= grid_size) continue;

                    int psf_idx = psf_r * grid_size + psf_c;
                    int pad_idx = (pad_start_y + rr) * P + (pad_start_x + cc);

                    lenslet_padded[pad_idx][0] = psf[psf_idx][0];
                    lenslet_padded[pad_idx][1] = psf[psf_idx][1];
                }
            }

            fftshift(lenslet_padded, lenslet_padded_shifted, P, P);
            fftw_execute(ifs_plan);
            fftshift(farfield_shifted, farfield, P, P);

            double scale = 1.0 / (P * P);
            for (int ii = 0; ii < P * P; ii++)
            {
                farfield[ii][0] *= scale;
                farfield[ii][1] *= scale;
            }

            if (pinhole_diam_ldl > 0.0)
            {
                double pin_diam_px = pinhole_diam_ldl * (P / (double)L);
                double l_cx = (P - 1) / 2.0;
                double l_cy = (P - 1) / 2.0;
                double pin_R = pin_diam_px / 2.0;

                for (int rr = 0; rr < P; rr++)
                {
                    double dy = rr - l_cy;
                    for (int cc = 0; cc < P; cc++)
                    {
                        double dx = cc - l_cx;
                        double dist = sqrt(dx * dx + dy * dy);
                        double transmission = 0.0;

                        if (dist < pin_R - 0.7072)
                        {
                            transmission = 1.0;
                        }
                        else if (dist > pin_R + 0.7072)
                        {
                            transmission = 0.0;
                        }
                        else
                        {
                            int count = 0;
                            int K = 20;
                            for (int ky = 0; ky < K; ky++)
                            {
                                double sy = rr - 0.5 + (ky + 0.5) / K - l_cy;
                                for (int kx = 0; kx < K; kx++)
                                {
                                    double sx = cc - 0.5 + (kx + 0.5) / K - l_cx;
                                    if (sx * sx + sy * sy <= pin_R * pin_R)
                                    {
                                        count++;
                                    }
                                }
                            }
                            transmission = (double)count / (K * K);
                        }
                        int far_idx = rr * P + cc;
                        farfield[far_idx][0] *= transmission;
                        farfield[far_idx][1] *= transmission;
                    }
                }
            }

            for (int rr = 0; rr < P; rr++)
            {
                int out_r = ly * P + rr;
                for (int cc = 0; cc < P; cc++)
                {
                    int out_c = lx * P + cc;
                    int far_idx = rr * P + cc;
                    int out_idx = out_r * out_size + out_c;

                    double r = farfield[far_idx][0];
                    double im = farfield[far_idx][1];
                    ifsraw_im_accum[out_idx] += r * r + im * im;
                }
            }
        }
    }
}

/* Run the Integral Field Spectrograph (IFS) simulation (single point source mode) */
int run_ifs_simulation(
    const fftw_complex *psf,
    int                 grid_size,
    double              diameter,
    double              lenslet_size_ld,
    int                 lenslet_count,
    int                 ifs_fft_size,
    const char         *ifs_complex_file,
    const char         *ifs_intensity_file,
    const char         *ifs_flux_file,
    double              pinhole_diam_ldl,
    const char         *ifs_mask_file,
    double              init_pupil_energy)
{
    double pixels_per_ld = (double)grid_size / diameter;
    double lenslet_size_px = lenslet_size_ld * pixels_per_ld;
    int L = (int)round(lenslet_size_px);

    if (L <= 0)
    {
        fprintf(stderr, "ERROR: [%s:%d] Calculated lenslet size is too small (%d pixels).\n",
                __func__, __LINE__, L);
        return 1;
    }
    if (L > ifs_fft_size)
    {
        fprintf(stderr, "ERROR: [%s:%d] Lenslet size (%d px) exceeds IFS FFT size (%d px).\n",
                __func__, __LINE__, L, ifs_fft_size);
        return 1;
    }

    int M = lenslet_count;
    int P = ifs_fft_size;
    int out_size = M * P;

    printf("IFS Simulation Configuration:\n");
    printf("  Lenslet size: %.3f lambda/D (%.2f pixels, rounded to %d px)\n",
           lenslet_size_ld, lenslet_size_px, L);
    printf("  Array grid size: %d x %d lenslets\n", M, M);
    printf("  Far-field FFT size: %d x %d (padded from %d px)\n", P, P, L);
    printf("  Arranged output grid size: %d x %d pixels\n", out_size, out_size);

    /* Load mask file if provided (applied to the tiled IFS output detector plane) */
    double *mask_data = NULL;
    if (ifs_mask_file != NULL)
    {
        int mask_rows = 0;
        int mask_cols = 0;
        printf("Loading IFS output mask from: %s\n", ifs_mask_file);
        mask_data = read_fits_2d(ifs_mask_file, &mask_rows, &mask_cols);
        if (mask_data == NULL)
        {
            fprintf(stderr, "ERROR: [%s:%d] Failed to load mask from %s.\n",
                    __func__, __LINE__, ifs_mask_file);
            return 1;
        }
        if (mask_rows != out_size || mask_cols != out_size)
        {
            fprintf(stderr, "ERROR: [%s:%d] Mask size (%d x %d) does not match "
                    "arranged IFS output size (%d x %d).\n",
                    __func__, __LINE__, mask_rows, mask_cols, out_size, out_size);
            free(mask_data);
            return 1;
        }
    }

    if (pinhole_diam_ldl > 0.0)
    {
        printf("  Pinhole mask enabled in IFS output: diameter %.3f lambda/dl\n", pinhole_diam_ldl);
    }

    /* Allocate output arrays */
    fftw_complex *ifs_complex = fftw_alloc_complex(out_size * out_size);
    double *ifs_intensity = calloc(out_size * out_size, sizeof(double));

    if (ifs_complex == NULL || ifs_intensity == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Out of memory for IFS arrays.\n", __func__, __LINE__);
        if (ifs_complex) fftw_free(ifs_complex);
        if (ifs_intensity) free(ifs_intensity);
        if (mask_data) free(mask_data);
        return 1;
    }

    /* Set up FFT buffers for a single lenslet */
    fftw_complex *lenslet_padded = fftw_alloc_complex(P * P);
    fftw_complex *lenslet_padded_shifted = fftw_alloc_complex(P * P);
    fftw_complex *farfield_shifted = fftw_alloc_complex(P * P);
    fftw_complex *farfield = fftw_alloc_complex(P * P);

    if (lenslet_padded == NULL || lenslet_padded_shifted == NULL ||
        farfield_shifted == NULL || farfield == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Out of memory for single-lenslet FFT buffers.\n",
                __func__, __LINE__);
        if (lenslet_padded) fftw_free(lenslet_padded);
        if (lenslet_padded_shifted) fftw_free(lenslet_padded_shifted);
        if (farfield_shifted) fftw_free(farfield_shifted);
        if (farfield) fftw_free(farfield);
        fftw_free(ifs_complex);
        free(ifs_intensity);
        if (mask_data) free(mask_data);
        return 1;
    }

    fftw_plan plan = fftw_plan_dft_2d(P, P, lenslet_padded_shifted, farfield_shifted,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
    if (!plan)
    {
        fprintf(stderr, "ERROR: [%s:%d] Failed to create FFTW plan for IFS.\n", __func__, __LINE__);
        fftw_free(lenslet_padded);
        fftw_free(lenslet_padded_shifted);
        fftw_free(farfield_shifted);
        fftw_free(farfield);
        fftw_free(ifs_complex);
        free(ifs_intensity);
        if (mask_data) free(mask_data);
        return 1;
    }

    /* Run core simulation loop */
    accumulate_ifs_source(psf, grid_size, L, M, P, lenslet_padded, lenslet_padded_shifted,
                          farfield_shifted, farfield, plan, pinhole_diam_ldl, ifs_intensity);

    /* Construct complex outputs and apply mask overlay if provided */
    for (int ii = 0; ii < out_size * out_size; ii++)
    {
        /* Reconstruct complex field from intensity
           (assuming zero phase for single source serialization) */
        ifs_complex[ii][0] = sqrt(ifs_intensity[ii]);
        ifs_complex[ii][1] = 0.0;
    }

    if (mask_data != NULL)
    {
        for (int ii = 0; ii < out_size * out_size; ii++)
        {
            double mask_val = mask_data[ii];
            ifs_complex[ii][0] *= mask_val;
            ifs_complex[ii][1] *= mask_val;
            ifs_intensity[ii] *= mask_val * mask_val;
        }
    }

    /* Calculate total throughput */
    double psf_energy = 0.0;
    for (int ii = 0; ii < grid_size * grid_size; ii++)
    {
        double r = psf[ii][0];
        double im = psf[ii][1];
        psf_energy += r * r + im * im;
    }

    double ifs_energy = 0.0;
    for (int ii = 0; ii < out_size * out_size; ii++)
    {
        ifs_energy += ifs_intensity[ii];
    }

    double throughput_ifs = 0.0;
    if (psf_energy > 0.0)
    {
        /* Multiply ifs_energy by P^2 to scale to physical energy unit of focal-plane PSF */
        throughput_ifs = (ifs_energy * P * P) / psf_energy;
    }

    double init_psf_energy = init_pupil_energy / (grid_size * grid_size);
    double throughput_system = 0.0;
    if (init_psf_energy > 0.0)
    {
        throughput_system = (ifs_energy * P * P) / init_psf_energy;
    }

    printf("  Throughput relative to post-coronagraph PSF: %.6e (%.3f%%)\n",
           throughput_ifs, throughput_ifs * 100.0);
    printf("  Total Coronagraph+IFS System Throughput:      "
           "%.6e (%.3f%% of initial pupil energy)\n",
           throughput_system, throughput_system * 100.0);

    /* Compute integrated flux per lenslet (M x M array) */
    double *ifs_flux = malloc(M * M * sizeof(double));
    if (ifs_flux == NULL)
    {
        fprintf(stderr, "ERROR: [%s:%d] Out of memory for ifs_flux array.\n", __func__, __LINE__);
    }
    else
    {
        double norm_scale = 0.0;
        if (init_pupil_energy > 0.0)
        {
            norm_scale = (double)(P * P) * (grid_size * grid_size) / init_pupil_energy;
        }

        for (int ly = 0; ly < M; ly++)
        {
            for (int lx = 0; lx < M; lx++)
            {
                double sum = 0.0;
                for (int rr = 0; rr < P; rr++)
                {
                    int out_r = ly * P + rr;
                    for (int cc = 0; cc < P; cc++)
                    {
                        int out_c = lx * P + cc;
                        sum += ifs_intensity[out_r * out_size + out_c];
                    }
                }
                ifs_flux[ly * M + lx] = sum * norm_scale;
            }
        }
        printf("Writing normalized IFS lenslet flux map to: %s\n", ifs_flux_file);
        write_fits_2d_wcs(ifs_flux_file, ifs_flux, M, M,
                          lenslet_size_ld, (M - 1) / 2.0 + 1.0);
        free(ifs_flux);
    }

    /* Write Output Files */
    printf("Writing IFS complex far-field array (amplitude, phase) to: %s\n", ifs_complex_file);
    write_fits_complex(ifs_complex_file, ifs_complex, out_size, out_size);

    printf("Writing IFS intensity array to: %s\n", ifs_intensity_file);
    write_fits_2d(ifs_intensity_file, ifs_intensity, out_size, out_size);

    /* Free buffers */
    fftw_destroy_plan(plan);
    fftw_free(lenslet_padded);
    fftw_free(lenslet_padded_shifted);
    fftw_free(farfield_shifted);
    fftw_free(farfield);
    fftw_free(ifs_complex);
    free(ifs_intensity);
    if (mask_data) free(mask_data);

    printf("IFS Simulation complete!\n");
    return 0;
}

/* Compute the focal plane complex amplitude of the matched reference mode */
void compute_mefs_reference_mode(
    const fftw_complex *pupil,
    int                 rows,
    int                 cols,
    double              diameter,
    double              coronagraph_param,
    double              mefs_x,
    double              mefs_y,
    fftw_complex       *mefs_mode)
{
    /* Copy base pupil */
    fftw_complex *curr_pupil = fftw_alloc_complex(rows * cols);
    for (int ii = 0; ii < rows * cols; ii++)
    {
        curr_pupil[ii][0] = pupil[ii][0];
        curr_pupil[ii][1] = pupil[ii][1];
    }

    /* Apply source tilt */
    if (mefs_x != 0.0 || mefs_y != 0.0)
    {
        double cx = (cols - 1) / 2.0;
        double cy = (rows - 1) / 2.0;
        for (int ii = 0; ii < rows; ii++)
        {
            double dy = ii - cy;
            for (int jj = 0; jj < cols; jj++)
            {
                double dx = jj - cx;
                int idx_p = ii * cols + jj;
                double phase_ramp = 2.0 * M_PI * (mefs_x * dx + mefs_y * dy) / diameter;
                double r = curr_pupil[idx_p][0];
                double im = curr_pupil[idx_p][1];
                curr_pupil[idx_p][0] = r * cos(phase_ramp) - im * sin(phase_ramp);
                curr_pupil[idx_p][1] = r * sin(phase_ramp) + im * cos(phase_ramp);
            }
        }
    }

    /* Apply coronagraph subtraction */
    double weights[6];
    get_coronagraph_weights(coronagraph_param, weights);

    int num_modes = 0;
    if (coronagraph_param > 0.0)
    {
        if (coronagraph_param <= 1.0) num_modes = 1;
        else if (coronagraph_param <= 2.0) num_modes = 3;
        else num_modes = 6;
    }

    if (num_modes > 0)
    {
        double G[36];
        double b_re[6];
        double b_im[6];
        memset(G, 0, sizeof(G));
        memset(b_re, 0, sizeof(b_re));
        memset(b_im, 0, sizeof(b_im));

        double cx = (cols - 1) / 2.0;
        double cy = (rows - 1) / 2.0;
        for (int ii = 0; ii < rows; ii++)
        {
            double dy = ii - cy;
            for (int jj = 0; jj < cols; jj++)
            {
                double dx = jj - cx;
                int idx = ii * cols + jj;
                double r = curr_pupil[idx][0];
                double im = curr_pupil[idx][1];
                double amp = sqrt(r * r + im * im);
                double u = 2.0 * dx / diameter;
                double v = 2.0 * dy / diameter;
                double val[6];
                val[0] = 1.0;
                val[1] = u;
                val[2] = v;
                val[3] = u * u;
                val[4] = u * v;
                val[5] = v * v;
                for (int m = 0; m < num_modes; m++)
                {
                    b_re[m] += r * val[m] * amp;
                    b_im[m] += im * val[m] * amp;
                    for (int k = 0; k < num_modes; k++)
                    {
                        G[m * num_modes + k] += val[m] * val[k] * amp * amp;
                    }
                }
            }
        }

        double G_copy[36];
        double c_re[6];
        double c_im[6];
        memcpy(c_re, b_re, sizeof(b_re));
        memcpy(c_im, b_im, sizeof(b_im));
        memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
        solve_linear_system(G_copy, c_re, num_modes);
        memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
        solve_linear_system(G_copy, c_im, num_modes);

        for (int ii = 0; ii < rows; ii++)
        {
            double dy = ii - cy;
            for (int jj = 0; jj < cols; jj++)
            {
                double dx = jj - cx;
                int idx_p = ii * cols + jj;
                double amp = sqrt(curr_pupil[idx_p][0] * curr_pupil[idx_p][0] +
                                  curr_pupil[idx_p][1] * curr_pupil[idx_p][1]);
                double u = 2.0 * dx / diameter;
                double v = 2.0 * dy / diameter;
                double val[6];
                val[0] = 1.0;
                val[1] = u;
                val[2] = v;
                val[3] = u * u;
                val[4] = u * v;
                val[5] = v * v;
                double sub_re = 0.0;
                double sub_im = 0.0;
                for (int m = 0; m < num_modes; m++)
                {
                    sub_re += weights[m] * c_re[m] * val[m] * amp;
                    sub_im += weights[m] * c_im[m] * val[m] * amp;
                }
                curr_pupil[idx_p][0] -= sub_re;
                curr_pupil[idx_p][1] -= sub_im;
            }
        }
    }

    /* Shift & FFT */
    fftw_complex *pupil_shifted = fftw_alloc_complex(rows * cols);
    fftw_complex *psf_shifted = fftw_alloc_complex(rows * cols);
    fftw_plan plan = fftw_plan_dft_2d(rows, cols, pupil_shifted, psf_shifted,
                                      FFTW_FORWARD, FFTW_ESTIMATE);

    fftshift(curr_pupil, pupil_shifted, rows, cols);
    fftw_execute(plan);
    fftshift(psf_shifted, mefs_mode, rows, cols);

    /* Scale complex amplitudes by 1 / (rows * cols) */
    double f_scale = 1.0 / (rows * cols);
    for (int ii = 0; ii < rows * cols; ii++)
    {
        mefs_mode[ii][0] *= f_scale;
        mefs_mode[ii][1] *= f_scale;
    }

    fftw_destroy_plan(plan);
    fftw_free(pupil_shifted);
    fftw_free(psf_shifted);
    fftw_free(curr_pupil);
}

int is_file_newer(
    const char *file_a,
    const char *file_b)
{
    struct stat stat_a;
    struct stat stat_b;
    if (stat(file_a, &stat_a) != 0)
    {
        return 0;
    }
    if (stat(file_b, &stat_b) != 0)
    {
        return 1;
    }
    return stat_a.st_mtime > stat_b.st_mtime;
}

void get_mefs_param_string(
    char   *buf,
    size_t  buf_len,
    double  coro_param,
    double  diameter,
    int     grid_size,
    double  mefs_x,
    double  mefs_y)
{
    snprintf(buf, buf_len,
             "# Params: coro_param=%.6f diameter=%.6f grid_size=%d mefs_x=%.6f mefs_y=%.6f",
             coro_param, diameter, grid_size, mefs_x, mefs_y);
}

int check_mefs_file_parameters(
    const char *filename,
    double      coro_param,
    double      diameter,
    int         grid_size,
    double      mefs_x,
    double      mefs_y)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL)
    {
        return 0;
    }

    char expected[256];
    get_mefs_param_string(expected, sizeof(expected), coro_param,
                          diameter, grid_size, mefs_x, mefs_y);

    char line[256];
    if (fgets(line, sizeof(line), f))
    {
        line[strcspn(line, "\r\n")] = '\0';
        if (strcmp(line, expected) == 0)
        {
            fclose(f);
            return 1;
        }
    }
    fclose(f);
    return 0;
}

void get_scene_param_string(
    char   *buf,
    size_t  buf_len,
    double  coro_param,
    double  diameter,
    int     grid_size,
    int     enable_ifs,
    double  lenslet_size,
    int     lenslet_count,
    int     ifs_fft_size,
    double  pinhole_diam)
{
    snprintf(buf, buf_len,
             "# Params: coro=%.6f diam=%.6f grid=%d ifs=%d lens_sz=%.6f "
             "lens_cnt=%d fft=%d pin_dia=%.6f",
             coro_param, diameter, grid_size, enable_ifs, lenslet_size,
             lenslet_count, ifs_fft_size, pinhole_diam);
}

int check_scene_file_parameters(
    const char *filename,
    double      coro_param,
    double      diameter,
    int         grid_size,
    int         enable_ifs,
    double      lenslet_size,
    int         lenslet_count,
    int         ifs_fft_size,
    double      pinhole_diam)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL)
    {
        return 0;
    }

    char expected[256];
    get_scene_param_string(expected, sizeof(expected), coro_param, diameter,
                           grid_size, enable_ifs, lenslet_size,
                           lenslet_count, ifs_fft_size, pinhole_diam);

    char line[256];
    if (fgets(line, sizeof(line), f))
    {
        line[strcspn(line, "\r\n")] = '\0';
        if (strcmp(line, expected) == 0)
        {
            fclose(f);
            return 1;
        }
    }
    fclose(f);
    return 0;
}

/* Load both uncoupled and coupled flux sums from a previously computed MEFS ASCII file */
int load_mefs_scene_data(
    const char *filename,
    double     *out_uncoupled,
    double     *out_coupled)
{
    FILE *f = fopen(filename, "r");
    if (f == NULL)
    {
        return 0;
    }

    double sum_uncoupled = 0.0;
    double sum_coupled = 0.0;
    char line[256];
    while (fgets(line, sizeof(line), f))
    {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
        {
            continue;
        }

        double x, y, flux, uncoupled, eff, coupled;
        if (sscanf(line, "%lf %lf %lf %lf %lf %lf", &x, &y, &flux, &uncoupled, &eff, &coupled) == 6)
        {
            sum_uncoupled += uncoupled;
            sum_coupled += coupled;
        }
    }
    fclose(f);
    *out_uncoupled = sum_uncoupled;
    *out_coupled = sum_coupled;
    return 1;
}

/* Light-weight scene simulation specifically optimized for MEFS coupling */
int run_mefs_scene_simulation(
    const char   *scene_path,
    const char   *name_suffix,
    fftw_complex *pupil,
    int           rows,
    int           cols,
    double        diameter,
    int           grid_size,
    double        coronagraph_param,
    double        mefs_x,
    double        mefs_y,
    double        init_pupil_energy,
    double       *out_total_uncoupled_flux,
    double       *out_total_coupled_flux)
{
    FILE *sf = fopen(scene_path, "r");
    if (sf == NULL)
    {
        fprintf(stderr, "ERROR: Could not open scene file %s\n", scene_path);
        return 1;
    }

    int scene_src_count = 0;
    char line[256];
    while (fgets(line, sizeof(line), sf))
    {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
        {
            continue;
        }

        double sx, sy, sflux;
        if (sscanf(line, "%lf %lf %lf", &sx, &sy, &sflux) == 3)
        {
            scene_src_count++;
        }
    }

    if (scene_src_count == 0)
    {
        fprintf(stderr, "ERROR: No valid sources found in scene file %s\n", scene_path);
        fclose(sf);
        return 1;
    }

    double *scene_src_x = malloc(scene_src_count * sizeof(double));
    double *scene_src_y = malloc(scene_src_count * sizeof(double));
    double *scene_src_flux = malloc(scene_src_count * sizeof(double));
    if (scene_src_x == NULL || scene_src_y == NULL || scene_src_flux == NULL)
    {
        fprintf(stderr, "ERROR: Out of memory for scene data.\n");
        if (scene_src_x) free(scene_src_x);
        if (scene_src_y) free(scene_src_y);
        if (scene_src_flux) free(scene_src_flux);
        fclose(sf);
        return 1;
    }

    rewind(sf);
    int idx = 0;
    double total_scene_flux = 0.0;
    while (fgets(line, sizeof(line), sf))
    {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r')
        {
            continue;
        }

        double sx, sy, sflux;
        if (sscanf(line, "%lf %lf %lf", &sx, &sy, &sflux) == 3)
        {
            scene_src_x[idx] = sx;
            scene_src_y[idx] = sy;
            scene_src_flux[idx] = sflux;
            total_scene_flux += sflux;
            idx++;
        }
    }
    fclose(sf);
    printf("Loaded %d sources from scene file %s (Total flux: %.6e)\n",
           scene_src_count, scene_path, total_scene_flux);

    /* Compute reference mode */
    fftw_complex *mefs_mode = fftw_alloc_complex(rows * cols);
    compute_mefs_reference_mode(pupil, rows, cols, diameter,
                                coronagraph_param, mefs_x, mefs_y, mefs_mode);

    double mefs_mode_norm = 0.0;
    for (int ii = 0; ii < rows * cols; ii++)
    {
        double r = mefs_mode[ii][0];
        double im = mefs_mode[ii][1];
        mefs_mode_norm += r * r + im * im;
    }

    int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp single
        num_threads = omp_get_num_threads();
    }
#endif

    fftw_complex **t_pupil_shifted = malloc(num_threads * sizeof(fftw_complex *));
    fftw_complex **t_psf_shifted = malloc(num_threads * sizeof(fftw_complex *));
    fftw_complex **t_psf = malloc(num_threads * sizeof(fftw_complex *));
    fftw_plan *t_plan = malloc(num_threads * sizeof(fftw_plan));

    for (int t = 0; t < num_threads; t++)
    {
        t_pupil_shifted[t] = fftw_alloc_complex(rows * cols);
        t_psf_shifted[t] = fftw_alloc_complex(rows * cols);
        t_psf[t] = fftw_alloc_complex(rows * cols);
        t_plan[t] = fftw_plan_dft_2d(rows, cols, t_pupil_shifted[t],
                                     t_psf_shifted[t], FFTW_FORWARD,
                                     FFTW_ESTIMATE);
    }

    double *coupled_fluxes = calloc(scene_src_count, sizeof(double));
    double *efficiencies = calloc(scene_src_count, sizeof(double));
    double *uncoupled_fluxes = calloc(scene_src_count, sizeof(double));

    int progress_counter = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int s = 0; s < scene_src_count; s++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif

        double sx = scene_src_x[s];
        double sy = scene_src_y[s];
        double sflux = scene_src_flux[s];

        fftw_complex *curr_pupil = fftw_alloc_complex(rows * cols);
        for (int ii = 0; ii < rows * cols; ii++)
        {
            curr_pupil[ii][0] = pupil[ii][0];
            curr_pupil[ii][1] = pupil[ii][1];
        }

        if (sx != 0.0 || sy != 0.0)
        {
            double cx = (cols - 1) / 2.0;
            double cy = (rows - 1) / 2.0;
            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx_p = ii * cols + jj;
                    double phase_ramp = 2.0 * M_PI *
                                        (sx * dx + sy * dy) / diameter;
                    double r = curr_pupil[idx_p][0];
                    double im = curr_pupil[idx_p][1];
                    curr_pupil[idx_p][0] = r * cos(phase_ramp) - im * sin(phase_ramp);
                    curr_pupil[idx_p][1] = r * sin(phase_ramp) + im * cos(phase_ramp);
                }
            }
        }

        double weights[6];
        get_coronagraph_weights(coronagraph_param, weights);

        int num_modes = 0;
        if (coronagraph_param > 0.0)
        {
            if (coronagraph_param <= 1.0) num_modes = 1;
            else if (coronagraph_param <= 2.0) num_modes = 3;
            else num_modes = 6;
        }

        if (num_modes > 0)
        {
            double G[36];
            double b_re[6];
            double b_im[6];
            memset(G, 0, sizeof(G));
            memset(b_re, 0, sizeof(b_re));
            memset(b_im, 0, sizeof(b_im));

            double cx = (cols - 1) / 2.0;
            double cy = (rows - 1) / 2.0;
            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx = ii * cols + jj;
                    double r = curr_pupil[idx][0];
                    double im = curr_pupil[idx][1];
                    double amp = sqrt(r * r + im * im);
                    double u = 2.0 * dx / diameter;
                    double v = 2.0 * dy / diameter;
                    double val[6];
                    val[0] = 1.0;
                    val[1] = u;
                    val[2] = v;
                    val[3] = u * u;
                    val[4] = u * v;
                    val[5] = v * v;
                    for (int m = 0; m < num_modes; m++)
                    {
                        b_re[m] += r * val[m] * amp;
                        b_im[m] += im * val[m] * amp;
                        for (int k = 0; k < num_modes; k++)
                        {
                            G[m * num_modes + k] += val[m] * val[k] * amp * amp;
                        }
                    }
                }
            }

            double G_copy[36];
            double c_re[6];
            double c_im[6];
            memcpy(c_re, b_re, sizeof(b_re));
            memcpy(c_im, b_im, sizeof(b_im));
            memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
            solve_linear_system(G_copy, c_re, num_modes);
            memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
            solve_linear_system(G_copy, c_im, num_modes);

            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx_p = ii * cols + jj;
                    double amp = sqrt(curr_pupil[idx_p][0] *
                                          curr_pupil[idx_p][0] +
                                      curr_pupil[idx_p][1] *
                                          curr_pupil[idx_p][1]);
                    double u = 2.0 * dx / diameter;
                    double v = 2.0 * dy / diameter;
                    double val[6];
                    val[0] = 1.0;
                    val[1] = u;
                    val[2] = v;
                    val[3] = u * u;
                    val[4] = u * v;
                    val[5] = v * v;
                    double sub_re = 0.0;
                    double sub_im = 0.0;
                    for (int m = 0; m < num_modes; m++)
                    {
                        sub_re += weights[m] * c_re[m] * val[m] * amp;
                        sub_im += weights[m] * c_im[m] * val[m] * amp;
                    }
                    curr_pupil[idx_p][0] -= sub_re;
                    curr_pupil[idx_p][1] -= sub_im;
                }
            }
        }

        fftshift(curr_pupil, t_pupil_shifted[tid], rows, cols);
        fftw_execute(t_plan[tid]);
        fftshift(t_psf_shifted[tid], t_psf[tid], rows, cols);

        double norm_scale = 1.0 / (rows * cols);
        for (int ii = 0; ii < rows * cols; ii++)
        {
            t_psf[tid][ii][0] *= norm_scale;
            t_psf[tid][ii][1] *= norm_scale;
        }

        double r_sum = 0.0;
        double i_sum = 0.0;
        for (int ii = 0; ii < rows * cols; ii++)
        {
            double ar = mefs_mode[ii][0];
            double ai = mefs_mode[ii][1];
            double br = t_psf[tid][ii][0];
            double bi = t_psf[tid][ii][1];
            r_sum += ar * br + ai * bi;
            i_sum += ar * bi - ai * br;
        }

        double eff = 0.0;
        if (mefs_mode_norm > 0.0)
        {
            double num = r_sum * r_sum + i_sum * i_sum;
            double den = mefs_mode_norm * mefs_mode_norm;
            eff = num / den;
        }

        double mefs_mode_norm_phys = 0.0;
        if (init_pupil_energy > 0.0)
        {
            mefs_mode_norm_phys = mefs_mode_norm * (grid_size * grid_size) /
                                  init_pupil_energy;
        }
        efficiencies[s] = eff;
        coupled_fluxes[s] = eff * sflux * mefs_mode_norm_phys;

        double post_coro_energy = init_pupil_energy;
        if (coronagraph_param > 0.0)
        {
            post_coro_energy = 0.0;
            for (int ii = 0; ii < rows * cols; ii++)
            {
                double r = curr_pupil[ii][0];
                double im = curr_pupil[ii][1];
                post_coro_energy += r * r + im * im;
            }
        }

        double uncoupled_flux = 0.0;
        if (init_pupil_energy > 0.0)
        {
            uncoupled_flux = sflux * (post_coro_energy / init_pupil_energy);
        }
        uncoupled_fluxes[s] = uncoupled_flux;

        fftw_free(curr_pupil);

        int my_progress;
#ifdef _OPENMP
#pragma omp atomic capture
#endif
        my_progress = ++progress_counter;

#ifdef _OPENMP
#pragma omp critical(progress_bar)
#endif
        {
            print_progress_bar(my_progress - 1, scene_src_count);
        }
    }
    printf("\n");

    char mefs_out[512];
    if (name_suffix != NULL)
    {
        snprintf(mefs_out, sizeof(mefs_out), "mefs.%s.txt", name_suffix);
    }
    else
    {
        snprintf(mefs_out, sizeof(mefs_out), "mefs.txt");
    }

    FILE *of = fopen(mefs_out, "w");
    if (of == NULL)
    {
        fprintf(stderr, "ERROR: Could not open output file %s\n", mefs_out);
    }
    else
    {
        char param_buf[256];
        get_mefs_param_string(param_buf, sizeof(param_buf), coronagraph_param,
                              diameter, grid_size, mefs_x, mefs_y);
        fprintf(of, "%s\n", param_buf);
        fprintf(of, "# MEFS Point Source Coupling Data\n");
        fprintf(of, "# Input Scene File: %s\n", scene_path);
        double total_uncoupled_flux = 0.0;
        double total_coupled_flux = 0.0;
        for (int s = 0; s < scene_src_count; s++)
        {
            total_uncoupled_flux += uncoupled_fluxes[s];
            total_coupled_flux += coupled_fluxes[s];
        }
        double avg_eff = 0.0;
        if (total_uncoupled_flux > 0.0)
        {
            avg_eff = total_coupled_flux / total_uncoupled_flux;
        }
        fprintf(of, "# Total Light in Focal Plane (before coupling): %.15e\n",
                total_uncoupled_flux);
        fprintf(of, "# Total Light Gathered:                         %.15e\n",
                total_coupled_flux);
        fprintf(of, "# Average Coupling Efficiency:                  %.15e\n",
                avg_eff);
        fprintf(of, "# Format:\n");
        fprintf(of, "# Col 1: X coordinate (lambda/D)\n");
        fprintf(of, "# Col 2: Y coordinate (lambda/D)\n");
        fprintf(of, "# Col 3: Input flux\n");
        fprintf(of, "# Col 4: Light in focal plane (before coupling)\n");
        fprintf(of, "# Col 5: Coupling efficiency (0.0 to 1.0)\n");
        fprintf(of, "# Col 6: Coupled flux\n");
        fprintf(of, "# ----------------------------------------------------------------------\n");

        for (int s = 0; s < scene_src_count; s++)
        {
            fprintf(of, "%10.4f %10.4f %16.8e %16.8e %16.8e %16.8e\n",
                    scene_src_x[s], scene_src_y[s], scene_src_flux[s],
                    uncoupled_fluxes[s], efficiencies[s], coupled_fluxes[s]);
        }
        fclose(of);
        printf("Successfully wrote MEFS coupling data to: %s\n", mefs_out);
    }

    double total_uncoupled_flux = 0.0;
    double total_coupled_flux = 0.0;
    for (int s = 0; s < scene_src_count; s++)
    {
        total_uncoupled_flux += uncoupled_fluxes[s];
        total_coupled_flux += coupled_fluxes[s];
    }
    *out_total_uncoupled_flux = total_uncoupled_flux;
    *out_total_coupled_flux = total_coupled_flux;

    free(efficiencies);
    free(coupled_fluxes);
    free(uncoupled_fluxes);
    fftw_free(mefs_mode);

    for (int t = 0; t < num_threads; t++)
    {
        fftw_destroy_plan(t_plan[t]);
        fftw_free(t_pupil_shifted[t]);
        fftw_free(t_psf_shifted[t]);
        fftw_free(t_psf[t]);
    }
    free(t_pupil_shifted);
    free(t_psf_shifted);
    free(t_psf);
    free(t_plan);

    free(scene_src_x);
    free(scene_src_y);
    free(scene_src_flux);

    return 0;
}

/* Helper function to run a point source scene simulation from file */
int run_scene_simulation_file(
    const char   *scene_path,
    const char   *name_suffix,
    fftw_complex *pupil,
    int           rows,
    int           cols,
    double        diameter,
    int           grid_size,
    int           enable_ifs,
    double        lenslet_size_ld,
    int           lenslet_count,
    int           ifs_fft_size,
    double        pinhole_diam_ldl,
    const char   *ifs_mask_file,
    double        coronagraph_param,
    double        init_pupil_energy,
    int           run_mefs_mode,
    double        mefs_x,
    double        mefs_y)
{
    char psf_out[512];
    char ifsraw_out[512];
    char ifs_out[512];
    char log_out[512];

    if (name_suffix != NULL)
    {
        snprintf(psf_out, sizeof(psf_out), "psf_im.%s.fits", name_suffix);
        snprintf(ifsraw_out, sizeof(ifsraw_out), "ifsraw_im.%s.fits", name_suffix);
        snprintf(ifs_out, sizeof(ifs_out), "ifs_im.%s.fits", name_suffix);
        snprintf(log_out, sizeof(log_out), "scene.%s.log", name_suffix);
    }
    else
    {
        snprintf(psf_out, sizeof(psf_out), "psf_im.fits");
        snprintf(ifsraw_out, sizeof(ifsraw_out), "ifsraw_im.fits");
        snprintf(ifs_out, sizeof(ifs_out), "ifs_im.fits");
        log_out[0] = '\0';
    }

    FILE *sf = fopen(scene_path, "r");
    if (sf == NULL)
    {
        fprintf(stderr, "ERROR: Could not open scene file %s\n", scene_path);
        return 1;
    }

    int scene_src_count = 0;
    char line[256];
    while (fgets(line, sizeof(line), sf))
    {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        double sx, sy, sflux;
        if (sscanf(line, "%lf %lf %lf", &sx, &sy, &sflux) == 3)
        {
            scene_src_count++;
        }
    }

    if (scene_src_count == 0)
    {
        fprintf(stderr, "ERROR: No valid sources found in scene file %s\n", scene_path);
        fclose(sf);
        return 1;
    }

    double *scene_src_x = malloc(scene_src_count * sizeof(double));
    double *scene_src_y = malloc(scene_src_count * sizeof(double));
    double *scene_src_flux = malloc(scene_src_count * sizeof(double));
    if (scene_src_x == NULL || scene_src_y == NULL || scene_src_flux == NULL)
    {
        fprintf(stderr, "ERROR: Out of memory for scene data.\n");
        if (scene_src_x) free(scene_src_x);
        if (scene_src_y) free(scene_src_y);
        if (scene_src_flux) free(scene_src_flux);
        fclose(sf);
        return 1;
    }

    rewind(sf);
    int idx = 0;
    double total_scene_flux = 0.0;
    while (fgets(line, sizeof(line), sf))
    {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        double sx, sy, sflux;
        if (sscanf(line, "%lf %lf %lf", &sx, &sy, &sflux) == 3)
        {
            scene_src_x[idx] = sx;
            scene_src_y[idx] = sy;
            scene_src_flux[idx] = sflux;
            total_scene_flux += sflux;
            idx++;
        }
    }
    fclose(sf);
    printf("Loaded %d sources from scene file %s (Total flux: %.6e)\n",
           scene_src_count, scene_path, total_scene_flux);

    /* Allocate accumulator arrays */
    double *psf_im_accum = calloc(rows * cols, sizeof(double));
    double *ifsraw_im_accum = NULL;
    int M = lenslet_count;
    int P = ifs_fft_size;
    int out_size = M * P;

    if (enable_ifs)
    {
        ifsraw_im_accum = calloc(out_size * out_size, sizeof(double));
        if (ifsraw_im_accum == NULL)
        {
            fprintf(stderr, "ERROR: Out of memory for ifsraw accumulator.\n");
            free(scene_src_x);
            free(scene_src_y);
            free(scene_src_flux);
            free(psf_im_accum);
            return 1;
        }
    }
    double pixels_per_ld = (double)grid_size / diameter;
    double lenslet_size_px = lenslet_size_ld * pixels_per_ld;
    int L = (int)round(lenslet_size_px);

    fftw_complex *mefs_mode = NULL;
    double mefs_mode_norm = 0.0;
    double *mefs_mode_psf_intensity = NULL;
    double *mefs_mode_ifsraw = NULL;
    double *mefs_mode_ifsflux = NULL;

    if (run_mefs_mode)
    {
        printf("Running MEFS projection mode (reference at x=%.3f, y=%.3f)...\n",
               mefs_x, mefs_y);
        mefs_mode = fftw_alloc_complex(rows * cols);
        compute_mefs_reference_mode(pupil, rows, cols, diameter,
                                    coronagraph_param, mefs_x, mefs_y, mefs_mode);

        mefs_mode_psf_intensity = calloc(rows * cols, sizeof(double));
        for (int ii = 0; ii < rows * cols; ii++)
        {
            double r = mefs_mode[ii][0];
            double im = mefs_mode[ii][1];
            mefs_mode_norm += r * r + im * im;
            mefs_mode_psf_intensity[ii] = r * r + im * im;
        }

        if (enable_ifs)
        {
            mefs_mode_ifsraw = calloc(out_size * out_size, sizeof(double));
            mefs_mode_ifsflux = calloc(M * M, sizeof(double));

            fftw_complex *t_lenslet_padded = fftw_alloc_complex(P * P);
            fftw_complex *t_lenslet_padded_shifted = fftw_alloc_complex(P * P);
            fftw_complex *t_farfield_shifted = fftw_alloc_complex(P * P);
            fftw_complex *t_farfield = fftw_alloc_complex(P * P);
            fftw_plan t_ifs_plan = fftw_plan_dft_2d(P, P, t_lenslet_padded_shifted,
                                                   t_farfield_shifted,
                                                   FFTW_FORWARD, FFTW_ESTIMATE);

            accumulate_ifs_source(mefs_mode, grid_size, L, M, P,
                                 t_lenslet_padded, t_lenslet_padded_shifted,
                                 t_farfield_shifted, t_farfield, t_ifs_plan,
                                 pinhole_diam_ldl, mefs_mode_ifsraw);

            fftw_destroy_plan(t_ifs_plan);
            fftw_free(t_lenslet_padded);
            fftw_free(t_lenslet_padded_shifted);
            fftw_free(t_farfield_shifted);
            fftw_free(t_farfield);

            double *mask_data = NULL;
            if (ifs_mask_file != NULL)
            {
                int mask_rows = 0;
                int mask_cols = 0;
                mask_data = read_fits_2d(ifs_mask_file, &mask_rows, &mask_cols);
                if (mask_data != NULL)
                {
                    for (int ii = 0; ii < out_size * out_size; ii++)
                    {
                        mefs_mode_ifsraw[ii] *= mask_data[ii] * mask_data[ii];
                    }
                    free(mask_data);
                }
            }

            double norm_scale = 0.0;
            if (init_pupil_energy > 0.0)
            {
                norm_scale = (double)(P * P) * (grid_size * grid_size) / init_pupil_energy;
            }

            for (int ly = 0; ly < M; ly++)
            {
                for (int lx = 0; lx < M; lx++)
                {
                    double sum = 0.0;
                    for (int rr = 0; rr < P; rr++)
                    {
                        int out_r = ly * P + rr;
                        for (int cc = 0; cc < P; cc++)
                        {
                            int out_c = lx * P + cc;
                            sum += mefs_mode_ifsraw[out_r * out_size + out_c];
                        }
                    }
                    mefs_mode_ifsflux[ly * M + lx] = sum * norm_scale;
                }
            }
        }
    }

    /* Set up thread-local environment for OpenMP parallelization */
    int num_threads = 1;
#ifdef _OPENMP
#pragma omp parallel
    {
#pragma omp single
        num_threads = omp_get_num_threads();
    }
#endif

    printf("Running scene simulation on %d thread(s)...\n", num_threads);
    /* Allocate thread-local accumulator arrays */
    double **thread_psf_accum = calloc(num_threads, sizeof(double *));
    double **thread_ifs_accum = calloc(num_threads, sizeof(double *));
    double *thread_coro_energy = calloc(num_threads, sizeof(double));
    double *thread_mefs_power = NULL;
    if (run_mefs_mode)
    {
        thread_mefs_power = calloc(num_threads, sizeof(double));
    }

    for (int t = 0; t < num_threads; t++)
    {
        thread_psf_accum[t] = calloc(rows * cols, sizeof(double));
        if (enable_ifs)
        {
            thread_ifs_accum[t] = calloc(out_size * out_size, sizeof(double));
        }
    }

    /* Allocate thread-local scratch arrays and FFTW plans */
    fftw_complex **t_pupil_shifted = malloc(num_threads * sizeof(fftw_complex *));
    fftw_complex **t_psf_shifted = malloc(num_threads * sizeof(fftw_complex *));
    fftw_complex **t_psf = malloc(num_threads * sizeof(fftw_complex *));
    fftw_plan *t_plan = malloc(num_threads * sizeof(fftw_plan));

    fftw_complex **t_lenslet_padded =
        malloc(num_threads * sizeof(fftw_complex *));
    fftw_complex **t_lenslet_padded_shifted =
        malloc(num_threads * sizeof(fftw_complex *));
    fftw_complex **t_farfield_shifted =
        malloc(num_threads * sizeof(fftw_complex *));
    fftw_complex **t_farfield =
        malloc(num_threads * sizeof(fftw_complex *));
    fftw_plan *t_ifs_plan = malloc(num_threads * sizeof(fftw_plan));

    for (int t = 0; t < num_threads; t++)
    {
        t_pupil_shifted[t] = fftw_alloc_complex(rows * cols);
        t_psf_shifted[t] = fftw_alloc_complex(rows * cols);
        t_psf[t] = fftw_alloc_complex(rows * cols);
        t_plan[t] = fftw_plan_dft_2d(rows, cols, t_pupil_shifted[t], t_psf_shifted[t],
                                     FFTW_FORWARD, FFTW_ESTIMATE);

        if (enable_ifs)
        {
            t_lenslet_padded[t] = fftw_alloc_complex(P * P);
            t_lenslet_padded_shifted[t] = fftw_alloc_complex(P * P);
            t_farfield_shifted[t] = fftw_alloc_complex(P * P);
            t_farfield[t] = fftw_alloc_complex(P * P);
            t_ifs_plan[t] = fftw_plan_dft_2d(P, P, t_lenslet_padded_shifted[t],
                                             t_farfield_shifted[t],
                                             FFTW_FORWARD, FFTW_ESTIMATE);
        }
    }

    int progress_counter = 0;

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (int s = 0; s < scene_src_count; s++)
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif

        double sx = scene_src_x[s];
        double sy = scene_src_y[s];
        double sflux = scene_src_flux[s];

        /* Copy base pupil */
        fftw_complex *curr_pupil = fftw_alloc_complex(rows * cols);
        for (int ii = 0; ii < rows * cols; ii++)
        {
            curr_pupil[ii][0] = pupil[ii][0];
            curr_pupil[ii][1] = pupil[ii][1];
        }

        /* Apply source tilt */
        if (sx != 0.0 || sy != 0.0)
        {
            double cx = (cols - 1) / 2.0;
            double cy = (rows - 1) / 2.0;
            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx_p = ii * cols + jj;
                    double phase_ramp = 2.0 * M_PI * (sx * dx + sy * dy) / diameter;
                    double r = curr_pupil[idx_p][0];
                    double im = curr_pupil[idx_p][1];
                    curr_pupil[idx_p][0] = r * cos(phase_ramp) - im * sin(phase_ramp);
                    curr_pupil[idx_p][1] = r * sin(phase_ramp) + im * cos(phase_ramp);
                }
            }
        }

        /* Apply idealized perfect coronagraph projection subtraction if enabled */
        double weights[6];
        get_coronagraph_weights(coronagraph_param, weights);

        int num_modes = 0;
        if (coronagraph_param > 0.0)
        {
            if (coronagraph_param <= 1.0) num_modes = 1;
            else if (coronagraph_param <= 2.0) num_modes = 3;
            else num_modes = 6;
        }

        if (num_modes > 0)
        {
            double G[36];
            double b_re[6];
            double b_im[6];

            memset(G, 0, sizeof(G));
            memset(b_re, 0, sizeof(b_re));
            memset(b_im, 0, sizeof(b_im));

            double cx = (cols - 1) / 2.0;
            double cy = (rows - 1) / 2.0;

            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx = ii * cols + jj;

                    double r = curr_pupil[idx][0];
                    double im = curr_pupil[idx][1];
                    double amp = sqrt(r * r + im * im);

                    double u = 2.0 * dx / diameter;
                    double v = 2.0 * dy / diameter;

                    double val[6];
                    val[0] = 1.0;
                    val[1] = u;
                    val[2] = v;
                    val[3] = u * u;
                    val[4] = u * v;
                    val[5] = v * v;

                    for (int m = 0; m < num_modes; m++)
                    {
                        b_re[m] += r * val[m] * amp;
                        b_im[m] += im * val[m] * amp;

                        for (int k = 0; k < num_modes; k++)
                        {
                            G[m * num_modes + k] += val[m] * val[k] * amp * amp;
                        }
                    }
                }
            }

            double G_copy[36];
            double c_re[6];
            double c_im[6];

            memcpy(c_re, b_re, sizeof(b_re));
            memcpy(c_im, b_im, sizeof(b_im));

            memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
            solve_linear_system(G_copy, c_re, num_modes);

            memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
            solve_linear_system(G_copy, c_im, num_modes);

            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx_p = ii * cols + jj;
                    double amp = sqrt(curr_pupil[idx_p][0] * curr_pupil[idx_p][0] +
                                      curr_pupil[idx_p][1] * curr_pupil[idx_p][1]);

                    double u = 2.0 * dx / diameter;
                    double v = 2.0 * dy / diameter;

                    double val[6];
                    val[0] = 1.0;
                    val[1] = u;
                    val[2] = v;
                    val[3] = u * u;
                    val[4] = u * v;
                    val[5] = v * v;

                    double sub_re = 0.0;
                    double sub_im = 0.0;
                    for (int m = 0; m < num_modes; m++)
                    {
                        sub_re += weights[m] * c_re[m] * val[m] * amp;
                        sub_im += weights[m] * c_im[m] * val[m] * amp;
                    }

                    curr_pupil[idx_p][0] -= sub_re;
                    curr_pupil[idx_p][1] -= sub_im;
                }
            }
        }

        /* Record post-coronagraph pupil energy for throughput calculation */
        double curr_coro_energy = 0.0;
        for (int ii = 0; ii < rows * cols; ii++)
        {
            double r = curr_pupil[ii][0];
            double im = curr_pupil[ii][1];
            curr_coro_energy += r * r + im * im;
        }
        thread_coro_energy[tid] += curr_coro_energy * sflux;

        /* Shift & FFT */
        fftshift(curr_pupil, t_pupil_shifted[tid], rows, cols);
        fftw_execute(t_plan[tid]);
        fftshift(t_psf_shifted[tid], t_psf[tid], rows, cols);

        /* Scale by 1/(N^2) and by sqrt(flux) */
        double f_scale = sqrt(sflux) / (rows * cols);
        for (int ii = 0; ii < rows * cols; ii++)
        {
            t_psf[tid][ii][0] *= f_scale;
            t_psf[tid][ii][1] *= f_scale;
        }

        if (run_mefs_mode)
        {
            double r_sum = 0.0;
            double i_sum = 0.0;
            for (int ii = 0; ii < rows * cols; ii++)
            {
                double ar = mefs_mode[ii][0];
                double ai = mefs_mode[ii][1];
                double br = t_psf[tid][ii][0];
                double bi = t_psf[tid][ii][1];
                r_sum += ar * br + ai * bi;
                i_sum += ar * bi - ai * br;
            }
            double coupled_power = 0.0;
            if (mefs_mode_norm > 0.0)
            {
                coupled_power = (r_sum * r_sum + i_sum * i_sum) / mefs_mode_norm;
            }
            thread_mefs_power[tid] += coupled_power;
        }

        /* Accumulate PSF intensity */
        for (int ii = 0; ii < rows * cols; ii++)
        {
            double r = t_psf[tid][ii][0];
            double im = t_psf[tid][ii][1];
            thread_psf_accum[tid][ii] += r * r + im * im;
        }

        /* Accumulate IFS if enabled */
        if (enable_ifs)
        {
            accumulate_ifs_source(
                t_psf[tid], grid_size, L, M, P,
                t_lenslet_padded[tid], t_lenslet_padded_shifted[tid],
                t_farfield_shifted[tid], t_farfield[tid], t_ifs_plan[tid],
                pinhole_diam_ldl, thread_ifs_accum[tid]);
        }

        fftw_free(curr_pupil);

        int my_progress;
#ifdef _OPENMP
#pragma omp atomic capture
#endif
        my_progress = ++progress_counter;

#ifdef _OPENMP
#pragma omp critical(progress_bar)
#endif
        {
            print_progress_bar(my_progress - 1, scene_src_count);
        }
    }
    printf("\n");

    /* Merge thread-local accumulators */
    double total_coro_pupil_energy = 0.0;
    for (int t = 0; t < num_threads; t++)
    {
        for (int ii = 0; ii < rows * cols; ii++)
        {
            psf_im_accum[ii] += thread_psf_accum[t][ii];
        }
        if (enable_ifs)
        {
            for (int ii = 0; ii < out_size * out_size; ii++)
            {
                ifsraw_im_accum[ii] += thread_ifs_accum[t][ii];
            }
        }
        total_coro_pupil_energy += thread_coro_energy[t];
    }

    double total_mefs_power = 0.0;
    if (run_mefs_mode)
    {
        for (int t = 0; t < num_threads; t++)
        {
            total_mefs_power += thread_mefs_power[t];
        }
        free(thread_mefs_power);
        printf("  MEFS Combined Coupled Power: %.6e\n", total_mefs_power);
    }

    /* Write PSF Intensity FITS File */
    printf("Writing scene PSF intensity to: %s\n", psf_out);
    write_fits_2d(psf_out, psf_im_accum, rows, cols);

    /* Report coronagraph throughput if enabled */
    double total_init_pupil_energy = init_pupil_energy * total_scene_flux;
    double coronagraph_throughput = 0.0;
    if (coronagraph_param > 0.0)
    {
        if (total_init_pupil_energy > 0.0)
        {
            coronagraph_throughput = total_coro_pupil_energy / total_init_pupil_energy;
        }
        printf("Coronagraph Throughput: %.6e (%.3f%% of initial pupil energy)\n",
               coronagraph_throughput, coronagraph_throughput * 100.0);
    }

    double scene_trans = 0.0;
    /* Finalize and write IFS files if enabled */
    if (enable_ifs)
    {
        /* Load mask overlay if provided */
        double *mask_data = NULL;
        if (ifs_mask_file != NULL)
        {
            int mask_rows = 0;
            int mask_cols = 0;
            mask_data = read_fits_2d(ifs_mask_file, &mask_rows, &mask_cols);
            if (mask_data != NULL)
            {
                for (int ii = 0; ii < out_size * out_size; ii++)
                {
                    ifsraw_im_accum[ii] *= mask_data[ii] * mask_data[ii];
                }
            }
        }

        /* Write raw IFS intensity */
        printf("Writing raw scene IFS intensity to: %s\n", ifsraw_out);
        write_fits_2d(ifsraw_out, ifsraw_im_accum, out_size, out_size);

        /* Calculate total throughput */
        double psf_energy = 0.0;
        for (int ii = 0; ii < grid_size * grid_size; ii++)
        {
            psf_energy += psf_im_accum[ii];
        }

        double ifs_energy = 0.0;
        for (int ii = 0; ii < out_size * out_size; ii++)
        {
            ifs_energy += ifsraw_im_accum[ii];
        }

        double throughput_ifs = 0.0;
        if (psf_energy > 0.0)
        {
            throughput_ifs = (ifs_energy * P * P) / psf_energy;
        }

        double init_psf_energy = total_init_pupil_energy / (grid_size * grid_size);
        double throughput_system = 0.0;
        if (init_psf_energy > 0.0)
        {
            throughput_system = (ifs_energy * P * P) / init_psf_energy;
        }

        printf("  Throughput relative to post-coronagraph PSF: %.6e (%.3f%%)\n",
               throughput_ifs, throughput_ifs * 100.0);
        printf("  Total Coronagraph+IFS System Throughput:      "
               "%.6e (%.3f%% of initial pupil energy)\n",
               throughput_system, throughput_system * 100.0);

        /* Compute normalized flux map */
        double *ifs_flux = malloc(M * M * sizeof(double));
        if (ifs_flux == NULL)
        {
            fprintf(stderr, "ERROR: Out of memory for ifs_flux array.\n");
        }
        else
        {
            double norm_scale = 0.0;
            if (init_pupil_energy > 0.0)
            {
                norm_scale = (double)(P * P) * (grid_size * grid_size) / init_pupil_energy;
            }

            double ifs_im_total_flux = 0.0;
            for (int ly = 0; ly < M; ly++)
            {
                for (int lx = 0; lx < M; lx++)
                {
                    double sum = 0.0;
                    for (int rr = 0; rr < P; rr++)
                    {
                        int out_r = ly * P + rr;
                        for (int cc = 0; cc < P; cc++)
                        {
                            int out_c = lx * P + cc;
                            sum += ifsraw_im_accum[out_r * out_size + out_c];
                        }
                    }
                    ifs_flux[ly * M + lx] = sum * norm_scale;
                    ifs_im_total_flux += ifs_flux[ly * M + lx];
                }
            }
            printf("Writing normalized scene IFS flux map to: %s\n", ifs_out);
            write_fits_2d_wcs(ifs_out, ifs_flux, M, M,
                              lenslet_size_ld, (M - 1) / 2.0 + 1.0);

            if (total_scene_flux > 0.0)
            {
                scene_trans = ifs_im_total_flux / total_scene_flux;
            }
            printf("  Total Flux in ifs_im.fits:                    %.6e\n",
                   ifs_im_total_flux);
            printf("  Total Input Scene Flux:                       %.6e\n",
                   total_scene_flux);
            printf("  Total Throughput (ifs_im / input):            %.6e (%.3f%%)\n",
                   scene_trans, scene_trans * 100.0);

            free(ifs_flux);
        }

        if (mask_data) free(mask_data);
    }

    /* Write summary log file if suffix is specified */
    if (name_suffix != NULL)
    {
        FILE *lf = fopen(log_out, "w");
        if (lf != NULL)
        {
            char param_buf[256];
            get_scene_param_string(param_buf, sizeof(param_buf), coronagraph_param,
                                  diameter, grid_size, enable_ifs, lenslet_size_ld,
                                  lenslet_count, ifs_fft_size, pinhole_diam_ldl);
            fprintf(lf, "%s\n", param_buf);
            fprintf(lf, "Scene Name: %s\n", name_suffix);
            fprintf(lf, "Input Scene File: %s\n", scene_path);
            fprintf(lf, "Total Input Scene Flux: %.15e\n", total_scene_flux);
            if (coronagraph_param > 0.0)
            {
                fprintf(lf, "Coronagraph Throughput: %.15e (%.6f%%)\n",
                        coronagraph_throughput, coronagraph_throughput * 100.0);
            }
            if (enable_ifs)
            {
                fprintf(lf, "IFS Throughput: %.15e (%.6f%%)\n",
                        scene_trans, scene_trans * 100.0);
            }
            if (run_mefs_mode)
            {
                fprintf(lf, "MEFS Combined Coupled Power: %.15e\n", total_mefs_power);
            }
            fclose(lf);
            printf("Saved summary log to: %s\n", log_out);
        }
    }

    /* Clean up thread-local arrays and plans */
    for (int t = 0; t < num_threads; t++)
    {
        free(thread_psf_accum[t]);
        fftw_destroy_plan(t_plan[t]);
        fftw_free(t_pupil_shifted[t]);
        fftw_free(t_psf_shifted[t]);
        fftw_free(t_psf[t]);

        if (enable_ifs)
        {
            free(thread_ifs_accum[t]);
            fftw_destroy_plan(t_ifs_plan[t]);
            fftw_free(t_lenslet_padded[t]);
            fftw_free(t_lenslet_padded_shifted[t]);
            fftw_free(t_farfield_shifted[t]);
            fftw_free(t_farfield[t]);
        }
    }
    free(thread_psf_accum);
    free(thread_ifs_accum);
    free(thread_coro_energy);

    free(t_pupil_shifted);
    free(t_psf_shifted);
    free(t_psf);
    free(t_plan);

    free(t_lenslet_padded);
    free(t_lenslet_padded_shifted);
    free(t_farfield_shifted);
    free(t_farfield);
    free(t_ifs_plan);

    /* Clean up shared allocations */
    free(psf_im_accum);
    if (enable_ifs) free(ifsraw_im_accum);
    free(scene_src_x);
    free(scene_src_y);
    free(scene_src_flux);

    if (run_mefs_mode)
    {
        fftw_free(mefs_mode);
        free(mefs_mode_psf_intensity);
        if (enable_ifs)
        {
            free(mefs_mode_ifsraw);
            free(mefs_mode_ifsflux);
        }
    }

    printf("Scene simulation complete!\n");
    return 0;
}

int main(
    int   argc,
    char **argv)
{
    init_colors();

    char *input_file = NULL;
    char *output_psf_file = "psf_complex.fits";
    char *output_pupil_file = "pupil.fits";
    char *output_intensity_file = "psf_intensity.fits";

    int enable_ifs = 0;
    double lenslet_size_ld = 0.5;
    int lenslet_count = 32;
    int ifs_fft_size = 32;
    char *ifs_complex_file = "ifs_complex.fits";
    char *ifs_intensity_file = "ifs_intensity.fits";
    char *ifs_flux_file = "ifs_flux.fits";
    double pinhole_diam_ldl = 2.0;
    char *ifs_mask_file = NULL;

    double source_x = 0.0;
    double source_y = 0.0;
    double coronagraph_param = 0.0;

    int grid_size = 512;
    double diameter = -1.0;

    char *scene_file = NULL;
    int run_scene_mode = 0;

    int run_mefs_mode = 0;
    double mefs_x = 1.0;
    double mefs_y = 0.0;

    static struct option long_options[] = {
        {"ifs", no_argument, 0, 's'},
        {"lenslet_size", required_argument, 0, 'l'},
        {"lenslet_count", required_argument, 0, 'm'},
        {"ifs_fft_size", required_argument, 0, 'P'},
        {"pinhole_diam", required_argument, 0, 'k'},
        {"ifs_mask", required_argument, 0, OPT_IFS_MASK},
        {"source_x", required_argument, 0, 'x'},
        {"source_y", required_argument, 0, 'y'},
        {"coronagraph", required_argument, 0, 'c'},
        {"ifs_complex", required_argument, 0, OPT_IFS_COMPLEX},
        {"ifs_intensity", required_argument, 0, OPT_IFS_INTENSITY},
        {"ifs_flux", required_argument, 0, OPT_IFS_FLUX},
        {"scene", optional_argument, 0, 'S'},
        {"mefs", no_argument, 0, 'M'},
        {"mefs_x", required_argument, 0, OPT_MEFS_X},
        {"mefs_y", required_argument, 0, OPT_MEFS_Y},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "i:o:p:t:n:d:sl:m:P:k:x:y:c:S::hM",
                             long_options, &option_index)) != -1)
    {
        switch (opt)
        {
            case 'i':
                input_file = optarg;
                break;
            case 'o':
                output_psf_file = optarg;
                break;
            case 'p':
                output_pupil_file = optarg;
                break;
            case 't':
                output_intensity_file = optarg;
                break;
            case 'n':
                grid_size = atoi(optarg);
                if (grid_size <= 0)
                {
                    fprintf(stderr, "ERROR: Grid size must be positive.\n");
                    return 1;
                }
                break;
            case 'd':
                diameter = atof(optarg);
                if (diameter <= 0)
                {
                    fprintf(stderr, "ERROR: Diameter must be positive.\n");
                    return 1;
                }
                break;
            case 's':
                enable_ifs = 1;
                break;
            case 'l':
                lenslet_size_ld = atof(optarg);
                if (lenslet_size_ld <= 0.0)
                {
                    fprintf(stderr, "ERROR: Lenslet size must be positive.\n");
                    return 1;
                }
                break;
            case 'm':
                lenslet_count = atoi(optarg);
                if (lenslet_count <= 0)
                {
                    fprintf(stderr, "ERROR: Lenslet count must be positive.\n");
                    return 1;
                }
                break;
            case 'P':
                ifs_fft_size = atoi(optarg);
                if (ifs_fft_size <= 0)
                {
                    fprintf(stderr, "ERROR: IFS FFT size must be positive.\n");
                    return 1;
                }
                break;
            case 'k':
                pinhole_diam_ldl = atof(optarg);
                if (pinhole_diam_ldl <= 0.0)
                {
                    fprintf(stderr, "ERROR: Pinhole diameter must be positive.\n");
                    return 1;
                }
                break;
            case 'x':
                source_x = atof(optarg);
                break;
            case 'y':
                source_y = atof(optarg);
                break;
            case 'c':
                coronagraph_param = atof(optarg);
                if (coronagraph_param < 0.0 || coronagraph_param > 3.0)
                {
                    fprintf(stderr,
                            "ERROR: Coronagraph param must be between 0.0 and 3.0.\n");
                    return 1;
                }
                break;
            case OPT_IFS_MASK:
                ifs_mask_file = optarg;
                break;
            case OPT_IFS_COMPLEX:
                ifs_complex_file = optarg;
                break;
            case OPT_IFS_INTENSITY:
                ifs_intensity_file = optarg;
                break;
            case OPT_IFS_FLUX:
                ifs_flux_file = optarg;
                break;
            case 'S':
                run_scene_mode = 1;
                if (optarg != NULL)
                {
                    scene_file = optarg;
                }
                else if (optind < argc && argv[optind][0] != '-')
                {
                    scene_file = argv[optind++];
                }
                else
                {
                    scene_file = NULL;
                }
                break;
            case 'M':
                run_mefs_mode = 1;
                break;
            case OPT_MEFS_X:
                mefs_x = atof(optarg);
                break;
            case OPT_MEFS_Y:
                mefs_y = atof(optarg);
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    int rows = grid_size;
    int cols = grid_size;
    fftw_complex *pupil = NULL;

    if (input_file != NULL)
    {
        printf("Reading pupil from FITS file: %s\n", input_file);
        pupil = read_fits_complex(input_file, &rows, &cols);
        if (pupil == NULL)
        {
            fprintf(stderr, "ERROR: Failed to read pupil from file.\n");
            return 1;
        }
        printf("Loaded pupil size: %d x %d\n", rows, cols);

        /* Estimate diameter if not provided */
        if (diameter < 0.0)
        {
            diameter = estimate_pupil_diameter(pupil, rows, cols);
            if (diameter < 0.0)
            {
                fprintf(stderr, "%sWARNING: Could not estimate pupil diameter. "
                        "Falling back to size/20.0%s\n", c_yellow, c_reset);
                diameter = grid_size / 20.0;
            }
            else
            {
                printf("Estimated pupil diameter from FITS file: %.2f pixels\n", diameter);
            }
        }
    }
    else
    {
        /* Set default diameter if not specified */
        if (diameter < 0)
        {
            diameter = grid_size / 20.0;
        }
        printf("Generating unobstructed circular pupil (grid: %d x %d, diameter: %.1f pix)\n",
               rows, cols, diameter);
        pupil = generate_circular_pupil(rows, cols, diameter);
        if (pupil == NULL)
        {
            return 1;
        }
    }

    /* Calculate initial pupil plane energy of unattenuated point source (for scaling) */
    double init_pupil_energy = 0.0;
    for (int ii = 0; ii < rows * cols; ii++)
    {
        double r = pupil[ii][0];
        double im = pupil[ii][1];
        init_pupil_energy += r * r + im * im;
    }

    if (run_scene_mode)
    {
        if (run_mefs_mode)
        {
            if (scene_file != NULL)
            {
                double total_uncoupled = 0.0;
                double total_coupled = 0.0;
                int status = run_mefs_scene_simulation(
                    scene_file, NULL, pupil, rows, cols, diameter,
                    grid_size, coronagraph_param, mefs_x, mefs_y,
                    init_pupil_energy, &total_uncoupled, &total_coupled);
                fftw_free(pupil);
                return status;
            }
            else
            {
                printf("Entering MEFS batch scene mode...\n");
                glob_t glob_results;
                int glob_ret = glob("scene.*.txt", 0, NULL, &glob_results);
                if (glob_ret != 0)
                {
                    printf("No scene files matching 'scene.*.txt' found.\n");
                    globfree(&glob_results);
                    fftw_free(pupil);
                    return 0;
                }

                char scene_names[100][256];
                double scene_uncoupled_fluxes[100];
                double scene_coupled_fluxes[100];
                int num_scenes = 0;

                for (size_t i = 0; i < glob_results.gl_pathc; i++)
                {
                    char *filename = glob_results.gl_pathv[i];
                    char *base = strrchr(filename, '/');
                    if (base == NULL) base = filename;
                    else base++;

                    if (strncmp(base, "scene.", 6) == 0)
                    {
                        char name[256];
                        strncpy(name, base + 6, sizeof(name) - 1);
                        name[sizeof(name) - 1] = '\0';
                        char *ext = strstr(name, ".txt");
                        if (ext != NULL)
                        {
                            *ext = '\0';
                        }

                        char mefs_out[512];
                        snprintf(mefs_out, sizeof(mefs_out), "mefs.%s.txt", name);

                        double uncoupled_flux = 0.0;
                        double coupled_flux = 0.0;
                        if (check_mefs_file_parameters(mefs_out, coronagraph_param,
                                                       diameter, grid_size, mefs_x, mefs_y) &&
                            !is_file_newer(filename, mefs_out))
                        {
                            printf("Skipping MEFS scene %s: cached parameters match. "
                                   "Loading existing coupling data...\n",
                                   name);
                            load_mefs_scene_data(mefs_out, &uncoupled_flux, &coupled_flux);
                        }
                        else
                        {
                            printf("Processing MEFS scene: %s (%s)...\n", name, filename);
                            run_mefs_scene_simulation(
                                filename, name, pupil, rows, cols, diameter,
                                grid_size, coronagraph_param, mefs_x, mefs_y,
                                init_pupil_energy, &uncoupled_flux, &coupled_flux);
                        }

                        strncpy(scene_names[num_scenes], name, 255);
                        scene_names[num_scenes][255] = '\0';
                        scene_uncoupled_fluxes[num_scenes] = uncoupled_flux;
                        scene_coupled_fluxes[num_scenes] = coupled_flux;
                        num_scenes++;
                    }
                }

                globfree(&glob_results);
                fftw_free(pupil);

                /* Generate MEFS Spectrograph Report */
                double planet_focal_light = 0.0;
                double planet_coupled_light = 0.0;
                double background_focal_light = 0.0;
                double background_coupled_light = 0.0;

                printf("\n========================================================\n");
                printf("MEFS Spectrograph Summary Report\n");
                printf("--------------------------------------------------------\n");

                for (int i = 0; i < num_scenes; i++)
                {
                    double u_flux = scene_uncoupled_fluxes[i];
                    double c_flux = scene_coupled_fluxes[i];
                    double eff = 0.0;
                    if (u_flux > 0.0)
                    {
                        eff = c_flux / u_flux;
                    }

                    printf("  %s scene:\n", scene_names[i]);
                    printf("    Light in Focal Plane (before coupling): %.6e\n", u_flux);
                    printf("    Coupled Light Gathered:                 %.6e\n", c_flux);
                    printf("    Average Coupling Efficiency:            %.6f (%.4f%%)\n",
                           eff, eff * 100.0);

                    if (strcmp(scene_names[i], "planet") == 0)
                    {
                        planet_focal_light = u_flux;
                        planet_coupled_light = c_flux;
                    }
                    else
                    {
                        background_focal_light += u_flux;
                        background_coupled_light += c_flux;
                    }
                }

                double total_coupled_light = planet_coupled_light + background_coupled_light;
                double snr = 0.0;
                if (total_coupled_light > 0.0)
                {
                    snr = planet_coupled_light / sqrt(total_coupled_light);
                }

                double planet_avg_eff = 0.0;
                if (planet_focal_light > 0.0)
                {
                    planet_avg_eff = planet_coupled_light / planet_focal_light;
                }

                double bg_avg_eff = 0.0;
                if (background_focal_light > 0.0)
                {
                    bg_avg_eff = background_coupled_light / background_focal_light;
                }

                printf("--------------------------------------------------------\n");
                printf("  Planet Light Gathered:      %.6e (avg eff: %.4f%%)\n",
                       planet_coupled_light, planet_avg_eff * 100.0);
                printf("  Total Background Light:     %.6e (avg eff: %.4f%%)\n",
                       background_coupled_light, bg_avg_eff * 100.0);
                printf("  Combined Spectrograph SNR:  %.6e\n", snr);
                printf("========================================================\n\n");

                /* Write to mefs.log */
                FILE *lf = fopen("mefs.log", "w");
                if (lf != NULL)
                {
                    fprintf(lf, "# MEFS Spectrograph Summary Report\n");
                    for (int i = 0; i < num_scenes; i++)
                    {
                        double u_flux = scene_uncoupled_fluxes[i];
                        double c_flux = scene_coupled_fluxes[i];
                        double eff = 0.0;
                        if (u_flux > 0.0)
                        {
                            eff = c_flux / u_flux;
                        }
                        fprintf(lf, "Scene: %s\n", scene_names[i]);
                        fprintf(lf, "  Focal Plane Light: %.15e\n", u_flux);
                        fprintf(lf, "  Coupled Light:     %.15e\n", c_flux);
                        fprintf(lf, "  Average Efficiency: %.15e\n", eff);
                    }
                    fprintf(lf, "Planet Light Gathered: %.15e\n", planet_coupled_light);
                    fprintf(lf, "Total Background Light: %.15e\n", background_coupled_light);
                    fprintf(lf, "Combined Spectrograph SNR: %.15e\n", snr);
                    fclose(lf);
                    printf("Successfully saved MEFS report to: mefs.log\n");
                }

                return 0;
            }
        }
        else
        {
            if (scene_file != NULL)
            {
                int status = run_scene_simulation_file(
                    scene_file, NULL, pupil, rows, cols, diameter,
                    grid_size, enable_ifs, lenslet_size_ld,
                    lenslet_count, ifs_fft_size, pinhole_diam_ldl,
                    ifs_mask_file, coronagraph_param, init_pupil_energy,
                    run_mefs_mode, mefs_x, mefs_y);
                fftw_free(pupil);
                return status;
            }
            else
            {
                printf("Entering batch scene simulation mode...\n");
                glob_t glob_results;
                int glob_ret = glob("scene.*.txt", 0, NULL, &glob_results);
                if (glob_ret != 0)
                {
                    printf("No scene files matching 'scene.*.txt' found.\n");
                    globfree(&glob_results);
                    fftw_free(pupil);
                    return 0;
                }

                /* Allocate master accumulators */
                double *master_psf_accum = calloc(rows * cols, sizeof(double));
                double *master_ifsraw_accum = NULL;
                double *master_ifs_flux_accum = NULL;
                int M = lenslet_count;
                int P = ifs_fft_size;
                int out_size = M * P;

                if (enable_ifs)
                {
                    master_ifsraw_accum = calloc(out_size * out_size, sizeof(double));
                    master_ifs_flux_accum = calloc(M * M, sizeof(double));
                }

                for (size_t i = 0; i < glob_results.gl_pathc; i++)
                {
                    char *filename = glob_results.gl_pathv[i];
                    char *base = strrchr(filename, '/');
                    if (base == NULL) base = filename;
                    else base++;

                    if (strncmp(base, "scene.", 6) == 0)
                    {
                        char name[256];
                        strncpy(name, base + 6, sizeof(name) - 1);
                        name[sizeof(name) - 1] = '\0';
                        char *ext = strstr(name, ".txt");
                        if (ext != NULL)
                        {
                            *ext = '\0';
                        }

                        char log_out[512];
                        snprintf(log_out, sizeof(log_out), "scene.%s.log", name);

                        char psf_out[512];
                        char ifsraw_out[512];
                        char ifs_out[512];
                        snprintf(psf_out, sizeof(psf_out), "psf_im.%s.fits", name);
                        snprintf(ifsraw_out, sizeof(ifsraw_out), "ifsraw_im.%s.fits", name);
                        snprintf(ifs_out, sizeof(ifs_out), "ifs_im.%s.fits", name);

                        if (check_scene_file_parameters(log_out, coronagraph_param,
                                                        diameter, grid_size, enable_ifs,
                                                        lenslet_size_ld, lenslet_count,
                                                        ifs_fft_size, pinhole_diam_ldl) &&
                            !is_file_newer(filename, log_out))
                        {
                            printf("Skipping scene %s: cached parameters match. "
                                   "Loading existing FITS images to sum...\n",
                                   name);

                            int fr = 0, fc = 0;
                            double *p_data = read_fits_2d(psf_out, &fr, &fc);
                            if (p_data != NULL)
                            {
                                if (fr == rows && fc == cols)
                                {
                                    for (int ii = 0; ii < rows * cols; ii++)
                                    {
                                        master_psf_accum[ii] += p_data[ii];
                                    }
                                }
                                free(p_data);
                            }

                            if (enable_ifs)
                            {
                                double *ir_data = read_fits_2d(ifsraw_out, &fr, &fc);
                                if (ir_data != NULL)
                                {
                                    if (fr == out_size && fc == out_size)
                                    {
                                        for (int ii = 0; ii < out_size * out_size; ii++)
                                        {
                                            master_ifsraw_accum[ii] += ir_data[ii];
                                        }
                                    }
                                    free(ir_data);
                                }

                                double *if_data = read_fits_2d(ifs_out, &fr, &fc);
                                if (if_data != NULL)
                                {
                                    if (fr == M && fc == M)
                                    {
                                        for (int ii = 0; ii < M * M; ii++)
                                        {
                                            master_ifs_flux_accum[ii] += if_data[ii];
                                        }
                                    }
                                    free(if_data);
                                }
                            }
                            continue;
                        }

                        printf("Processing batch scene: %s (%s)...\n", name, filename);
                        run_scene_simulation_file(
                            filename, name, pupil, rows, cols, diameter,
                            grid_size, enable_ifs, lenslet_size_ld,
                            lenslet_count, ifs_fft_size, pinhole_diam_ldl,
                            ifs_mask_file, coronagraph_param, init_pupil_energy,
                            run_mefs_mode, mefs_x, mefs_y);

                        /* Read and accumulate from the newly written files */
                        int fr = 0, fc = 0;
                        double *p_data = read_fits_2d(psf_out, &fr, &fc);
                        if (p_data != NULL)
                        {
                            if (fr == rows && fc == cols)
                            {
                                for (int ii = 0; ii < rows * cols; ii++)
                                {
                                    master_psf_accum[ii] += p_data[ii];
                                }
                            }
                            free(p_data);
                        }

                        if (enable_ifs)
                        {
                            double *ir_data = read_fits_2d(ifsraw_out, &fr, &fc);
                            if (ir_data != NULL)
                            {
                                if (fr == out_size && fc == out_size)
                                {
                                    for (int ii = 0; ii < out_size * out_size; ii++)
                                    {
                                        master_ifsraw_accum[ii] += ir_data[ii];
                                    }
                                }
                                free(ir_data);
                            }

                            double *if_data = read_fits_2d(ifs_out, &fr, &fc);
                            if (if_data != NULL)
                            {
                                if (fr == M && fc == M)
                                {
                                    for (int ii = 0; ii < M * M; ii++)
                                    {
                                        master_ifs_flux_accum[ii] += if_data[ii];
                                    }
                                }
                                free(if_data);
                            }
                        }
                    }
                }

                /* Write the combined sum images */
                printf("Writing combined scene PSF sum to: psf_im.scenesum.fits\n");
                write_fits_2d("psf_im.scenesum.fits", master_psf_accum, rows, cols);
                free(master_psf_accum);

                if (enable_ifs)
                {
                    printf("Writing combined raw scene IFS sum to: ifsraw_im.scenesum.fits\n");
                    write_fits_2d("ifsraw_im.scenesum.fits", master_ifsraw_accum,
                                  out_size, out_size);
                    free(master_ifsraw_accum);

                    printf("Writing combined normalized scene IFS flux sum to: "
                           "ifs_im.scenesum.fits\n");
                    write_fits_2d_wcs("ifs_im.scenesum.fits", master_ifs_flux_accum, M, M,
                                      lenslet_size_ld, (M - 1) / 2.0 + 1.0);
                    free(master_ifs_flux_accum);
                }

                globfree(&glob_results);
                fftw_free(pupil);
                printf("Batch scene simulation complete.\n");
                return 0;
            }
        }
    }

    if (run_mefs_mode)
    {
        /* Single source MEFS mode - no FITS files created */
        fftw_complex *untilted_pupil_mefs = fftw_alloc_complex(rows * cols);
        if (untilted_pupil_mefs != NULL)
        {
            for (int ii = 0; ii < rows * cols; ii++)
            {
                untilted_pupil_mefs[ii][0] = pupil[ii][0];
                untilted_pupil_mefs[ii][1] = pupil[ii][1];
            }
        }

        if (source_x != 0.0 || source_y != 0.0)
        {
            printf("Applying off-axis source tilt: x = %.3f, y = %.3f lambda/D\n",
                   source_x, source_y);
            double cx = (cols - 1) / 2.0;
            double cy = (rows - 1) / 2.0;
            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx = ii * cols + jj;
                    double phase_ramp = 2.0 * M_PI *
                                        (source_x * dx + source_y * dy) / diameter;
                    double r = pupil[idx][0];
                    double im = pupil[idx][1];
                    pupil[idx][0] = r * cos(phase_ramp) - im * sin(phase_ramp);
                    pupil[idx][1] = r * sin(phase_ramp) + im * cos(phase_ramp);
                }
            }
        }

        double weights[6];
        get_coronagraph_weights(coronagraph_param, weights);

        int num_modes = 0;
        if (coronagraph_param > 0.0)
        {
            if (coronagraph_param <= 1.0) num_modes = 1;
            else if (coronagraph_param <= 2.0) num_modes = 3;
            else num_modes = 6;
        }

        if (num_modes > 0)
        {

            double G[36];
            double b_re[6];
            double b_im[6];
            memset(G, 0, sizeof(G));
            memset(b_re, 0, sizeof(b_re));
            memset(b_im, 0, sizeof(b_im));

            double cx = (cols - 1) / 2.0;
            double cy = (rows - 1) / 2.0;
            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx = ii * cols + jj;
                    double r = pupil[idx][0];
                    double im = pupil[idx][1];
                    double amp = sqrt(r * r + im * im);
                    double u = 2.0 * dx / diameter;
                    double v = 2.0 * dy / diameter;
                    double val[6];
                    val[0] = 1.0;
                    val[1] = u;
                    val[2] = v;
                    val[3] = u * u;
                    val[4] = u * v;
                    val[5] = v * v;
                    for (int m = 0; m < num_modes; m++)
                    {
                        b_re[m] += r * val[m] * amp;
                        b_im[m] += im * val[m] * amp;
                        for (int k = 0; k < num_modes; k++)
                        {
                            G[m * num_modes + k] += val[m] * val[k] * amp * amp;
                        }
                    }
                }
            }

            double G_copy[36];
            double c_re[6];
            double c_im[6];
            memcpy(c_re, b_re, sizeof(b_re));
            memcpy(c_im, b_im, sizeof(b_im));
            memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
            solve_linear_system(G_copy, c_re, num_modes);
            memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
            solve_linear_system(G_copy, c_im, num_modes);

            for (int ii = 0; ii < rows; ii++)
            {
                double dy = ii - cy;
                for (int jj = 0; jj < cols; jj++)
                {
                    double dx = jj - cx;
                    int idx_p = ii * cols + jj;
                    double amp = sqrt(pupil[idx_p][0] * pupil[idx_p][0] +
                                      pupil[idx_p][1] * pupil[idx_p][1]);
                    double u = 2.0 * dx / diameter;
                    double v = 2.0 * dy / diameter;
                    double val[6];
                    val[0] = 1.0;
                    val[1] = u;
                    val[2] = v;
                    val[3] = u * u;
                    val[4] = u * v;
                    val[5] = v * v;
                    double sub_re = 0.0;
                    double sub_im = 0.0;
                    for (int m = 0; m < num_modes; m++)
                    {
                        sub_re += weights[m] * c_re[m] * val[m] * amp;
                        sub_im += weights[m] * c_im[m] * val[m] * amp;
                    }
                    pupil[idx_p][0] -= sub_re;
                    pupil[idx_p][1] -= sub_im;
                }
            }
        }

        fftw_complex *pupil_shifted = fftw_alloc_complex(rows * cols);
        fftw_complex *psf_shifted = fftw_alloc_complex(rows * cols);
        fftw_complex *psf = fftw_alloc_complex(rows * cols);
        fftw_plan plan = fftw_plan_dft_2d(rows, cols, pupil_shifted, psf_shifted,
                                          FFTW_FORWARD, FFTW_ESTIMATE);

        fftshift(pupil, pupil_shifted, rows, cols);
        fftw_execute(plan);
        fftshift(psf_shifted, psf, rows, cols);

        double scale = 1.0 / (rows * cols);
        for (int ii = 0; ii < rows * cols; ii++)
        {
            psf[ii][0] *= scale;
            psf[ii][1] *= scale;
        }

        printf("Running MEFS projection mode (reference at x=%.3f, y=%.3f)...\n",
               mefs_x, mefs_y);
        fftw_complex *mefs_mode = fftw_alloc_complex(rows * cols);
        compute_mefs_reference_mode(untilted_pupil_mefs, rows, cols, diameter,
                                    coronagraph_param, mefs_x, mefs_y, mefs_mode);

        double mefs_mode_norm = 0.0;
        double r_sum = 0.0;
        double i_sum = 0.0;
        for (int ii = 0; ii < rows * cols; ii++)
        {
            double ar = mefs_mode[ii][0];
            double ai = mefs_mode[ii][1];
            double br = psf[ii][0];
            double bi = psf[ii][1];
            r_sum += ar * br + ai * bi;
            i_sum += ar * bi - ai * br;
            mefs_mode_norm += ar * ar + ai * ai;
        }

        double c_re = 0.0;
        double c_im = 0.0;
        if (mefs_mode_norm > 0.0)
        {
            c_re = r_sum / mefs_mode_norm;
            c_im = i_sum / mefs_mode_norm;
        }
        double coupled_efficiency = c_re * c_re + c_im * c_im;
        double mefs_mode_norm_phys = 0.0;
        if (init_pupil_energy > 0.0)
        {
            mefs_mode_norm_phys = mefs_mode_norm * (rows * cols) / init_pupil_energy;
        }
        double coupled_flux = coupled_efficiency * mefs_mode_norm_phys;
        printf("  Coupling efficiency:    %.6e (%.4f%%)\n",
               coupled_efficiency, coupled_efficiency * 100.0);
        printf("  Coupled flux (flux=1):  %.6e\n", coupled_flux);

        fftw_free(mefs_mode);
        fftw_destroy_plan(plan);
        fftw_free(pupil);
        fftw_free(pupil_shifted);
        fftw_free(psf_shifted);
        fftw_free(psf);
        if (untilted_pupil_mefs != NULL) fftw_free(untilted_pupil_mefs);
        fftw_cleanup();
        printf("Done!\n");
        return 0;
    }

    /* --- SINGLE POINT SOURCE MODE (ORIGINAL FLOW) --- */

    fftw_complex *untilted_pupil = fftw_alloc_complex(rows * cols);
    if (untilted_pupil != NULL)
    {
        for (int ii = 0; ii < rows * cols; ii++)
        {
            untilted_pupil[ii][0] = pupil[ii][0];
            untilted_pupil[ii][1] = pupil[ii][1];
        }
    }

    /* Apply off-axis source wavefront tilt if specified */
    if (source_x != 0.0 || source_y != 0.0)
    {
        printf("Applying off-axis source tilt: x = %.3f, y = %.3f lambda/D\n", source_x, source_y);
        double cx = (cols - 1) / 2.0;
        double cy = (rows - 1) / 2.0;
        for (int ii = 0; ii < rows; ii++)
        {
            double dy = ii - cy;
            for (int jj = 0; jj < cols; jj++)
            {
                double dx = jj - cx;
                int idx = ii * cols + jj;
                double phase_ramp = 2.0 * M_PI * (source_x * dx + source_y * dy) / diameter;
                double r = pupil[idx][0];
                double im = pupil[idx][1];
                pupil[idx][0] = r * cos(phase_ramp) - im * sin(phase_ramp);
                pupil[idx][1] = r * sin(phase_ramp) + im * cos(phase_ramp);
            }
        }
    }

    /* Apply idealized perfect coronagraph subtraction if enabled */
    double weights[6];
    get_coronagraph_weights(coronagraph_param, weights);

    int num_modes = 0;
    if (coronagraph_param > 0.0)
    {
        printf("Applying idealized perfect coronagraph projection subtraction "
               "(Param %.6f)...\n", coronagraph_param);
        if (coronagraph_param <= 1.0) num_modes = 1;
        else if (coronagraph_param <= 2.0) num_modes = 3;
        else num_modes = 6;
    }

    if (num_modes > 0)
    {

        double G[36];
        double b_re[6];
        double b_im[6];

        memset(G, 0, sizeof(G));
        memset(b_re, 0, sizeof(b_re));
        memset(b_im, 0, sizeof(b_im));

        double cx = (cols - 1) / 2.0;
        double cy = (rows - 1) / 2.0;

        /* Calculate projection matrix G and vector b */
        for (int ii = 0; ii < rows; ii++)
        {
            double dy = ii - cy;
            for (int jj = 0; jj < cols; jj++)
            {
                double dx = jj - cx;
                int idx = ii * cols + jj;

                double r = pupil[idx][0];
                double im = pupil[idx][1];
                double amp = sqrt(r * r + im * im);

                double u = 2.0 * dx / diameter;
                double v = 2.0 * dy / diameter;

                double val[6];
                val[0] = 1.0;
                val[1] = u;
                val[2] = v;
                val[3] = u * u;
                val[4] = u * v;
                val[5] = v * v;

                for (int m = 0; m < num_modes; m++)
                {
                    b_re[m] += r * val[m] * amp;
                    b_im[m] += im * val[m] * amp;

                    for (int k = 0; k < num_modes; k++)
                    {
                        G[m * num_modes + k] += val[m] * val[k] * amp * amp;
                    }
                }
            }
        }

        /* Solve G * c = b for real and imaginary parts */
        double G_copy[36];
        double c_re[6];
        double c_im[6];

        memcpy(c_re, b_re, sizeof(b_re));
        memcpy(c_im, b_im, sizeof(b_im));

        memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
        if (solve_linear_system(G_copy, c_re, num_modes) != 0)
        {
            fprintf(stderr, "ERROR: Failed to solve real part of coronagraph coefficients.\n");
            return 1;
        }

        memcpy(G_copy, G, num_modes * num_modes * sizeof(double));
        if (solve_linear_system(G_copy, c_im, num_modes) != 0)
        {
            fprintf(stderr, "ERROR: Failed to solve imaginary part of coronagraph coefficients.\n");
            return 1;
        }

        printf("  Coronagraph subtraction coefficients:\n");
        for (int m = 0; m < num_modes; m++)
        {
            printf("    Mode %d: re = %.6f, im = %.6f\n", m, c_re[m], c_im[m]);
        }

        /* Subtract projection from pupil */
        for (int ii = 0; ii < rows; ii++)
        {
            double dy = ii - cy;
            for (int jj = 0; jj < cols; jj++)
            {
                double dx = jj - cx;
                int idx = ii * cols + jj;
                double amp = sqrt(pupil[idx][0] * pupil[idx][0] + pupil[idx][1] * pupil[idx][1]);

                double u = 2.0 * dx / diameter;
                double v = 2.0 * dy / diameter;

                double val[6];
                val[0] = 1.0;
                val[1] = u;
                val[2] = v;
                val[3] = u * u;
                val[4] = u * v;
                val[5] = v * v;

                double sub_re = 0.0;
                double sub_im = 0.0;
                for (int m = 0; m < num_modes; m++)
                {
                    sub_re += weights[m] * c_re[m] * val[m] * amp;
                    sub_im += weights[m] * c_im[m] * val[m] * amp;
                }

                pupil[idx][0] -= sub_re;
                pupil[idx][1] -= sub_im;
            }
        }

        /* Calculate post-coronagraph pupil plane energy */
        double coro_pupil_energy = 0.0;
        for (int ii = 0; ii < rows * cols; ii++)
        {
            double r = pupil[ii][0];
            double im = pupil[ii][1];
            coro_pupil_energy += r * r + im * im;
        }

        double coronagraph_throughput = 0.0;
        if (init_pupil_energy > 0.0)
        {
            coronagraph_throughput = coro_pupil_energy / init_pupil_energy;
        }
        printf("Coronagraph Throughput: %.6e (%.3f%% of initial pupil energy)\n",
               coronagraph_throughput, coronagraph_throughput * 100.0);
    }

    /* Saving pupil wavefront (with tilt and coronagraph subtraction if applied) */
    printf("Saving pupil to FITS file: %s\n", output_pupil_file);
    write_fits_complex(output_pupil_file, pupil, rows, cols);

    /* Allocate arrays for shifting and FFT */
    fftw_complex *pupil_shifted = fftw_alloc_complex(rows * cols);
    fftw_complex *psf_shifted = fftw_alloc_complex(rows * cols);
    fftw_complex *psf = fftw_alloc_complex(rows * cols);

    if (pupil_shifted == NULL || psf_shifted == NULL || psf == NULL)
    {
        fprintf(stderr, "ERROR: Failed to allocate arrays for FFT operations.\n");
        if (pupil != NULL) fftw_free(pupil);
        if (pupil_shifted != NULL) fftw_free(pupil_shifted);
        if (psf_shifted != NULL) fftw_free(psf_shifted);
        if (psf != NULL) fftw_free(psf);
        return 1;
    }

    /* Create FFTW Plan */
    fftw_plan plan = fftw_plan_dft_2d(rows, cols, pupil_shifted, psf_shifted,
                                      FFTW_FORWARD, FFTW_ESTIMATE);
    if (!plan)
    {
        fprintf(stderr, "ERROR: Failed to create FFTW plan.\n");
        fftw_free(pupil);
        fftw_free(pupil_shifted);
        fftw_free(psf_shifted);
        fftw_free(psf);
        return 1;
    }

    /* 1. Shift the input pupil from centered space to FFT zero-frequency corner space */
    fftshift(pupil, pupil_shifted, rows, cols);

    /* 2. Perform DFT */
    printf("Computing FFT...\n");
    fftw_execute(plan);

    /* 3. Shift zero-frequency back to the center of the grid */
    fftshift(psf_shifted, psf, rows, cols);

    /* 4. Normalize the PSF. Scale complex amplitudes by 1 / (rows * cols). */
    double scale = 1.0 / (rows * cols);
    for (int ii = 0; ii < rows * cols; ii++)
    {
        psf[ii][0] *= scale;
        psf[ii][1] *= scale;
    }



    /* 5. Compute Intensity PSF */
    double *intensity = malloc(rows * cols * sizeof(double));
    if (intensity != NULL)
    {
        for (int ii = 0; ii < rows * cols; ii++)
        {
            double r = psf[ii][0];
            double im = psf[ii][1];
            intensity[ii] = r * r + im * im;
        }
    }
    else
    {
        fprintf(stderr, "WARNING: Failed to allocate memory for intensity PSF.\n");
    }

    /* Write Output Files */
    printf("Writing complex amplitude PSF (amplitude, phase) to: %s\n", output_psf_file);
    write_fits_complex(output_psf_file, psf, rows, cols);

    if (intensity != NULL)
    {
        printf("Writing intensity PSF to: %s\n", output_intensity_file);
        write_fits_2d(output_intensity_file, intensity, rows, cols);
        free(intensity);
    }

    /* 6. Run IFS simulation if enabled */
    if (enable_ifs)
    {
        int ifs_status = run_ifs_simulation(psf, grid_size, diameter, lenslet_size_ld,
                                            lenslet_count, ifs_fft_size, ifs_complex_file,
                                            ifs_intensity_file, ifs_flux_file, pinhole_diam_ldl,
                                            ifs_mask_file, init_pupil_energy);
        if (ifs_status != 0)
        {
            fprintf(stderr, "ERROR: IFS simulation failed.\n");
            fftw_destroy_plan(plan);
            fftw_free(pupil);
            fftw_free(pupil_shifted);
            fftw_free(psf_shifted);
            fftw_free(psf);
            fftw_free(untilted_pupil);
            return 1;
        }
    }

    /* Clean up */
    fftw_destroy_plan(plan);
    fftw_free(pupil);
    fftw_free(pupil_shifted);
    fftw_free(psf_shifted);
    fftw_free(psf);
    fftw_free(untilted_pupil);
    fftw_cleanup();

    printf("Done!\n");
    return 0;
}
