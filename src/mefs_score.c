#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <getopt.h>
#include <fitsio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

const char *c_reset   = "";
const char *c_bold    = "";
const char *c_cyan    = "";
const char *c_b_green = "";
const char *c_green   = "";
const char *c_magenta = "";
const char *c_yellow  = "";
const char *c_b_red   = "";
const char *c_grey    = "";

typedef struct
{
    int    row;
    int    col;
    double x_ld;
    double y_ld;
    double val_planet;
    double val_scenesum;
    double val_diff;
} PixelScore;

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

void print_help(
    const char *prog_name)
{
    printf("%sNAME%s\n", c_cyan, c_reset);
    printf("  %s%s%s - Grade planet detection contrast metric from IFS outputs.\n\n",
           c_b_green, prog_name, c_reset);

    printf("%sUSAGE%s\n", c_cyan, c_reset);
    printf("  %s%s%s [OPTIONS]\n\n", c_b_green, prog_name, c_reset);

    printf("%sDESCRIPTION%s\n", c_cyan, c_reset);
    printf("  For each pixel in the planet image, measures total, planet, and background\n");
    printf("  flux, and outputs a sorted list by planet intensity. The output ASCII\n");
    printf("  file contains cumulative integrations to evaluate detection metrics.\n");
    printf("  The program automatically extracts the pixel sampling scale (CDELT1) and\n");
    printf("  origin (CRPIX1, CRPIX2) from the planet image FITS WCS header cards.\n\n");

    printf("%sOPTIONS%s\n", c_cyan, c_reset);
    printf("  %s-p, --planet%s %s<file>%s        "
           "Planet flux map FITS file (default: ifs_im.planet.fits).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-s, --scenesum%s %s<file>%s      "
           "Combined scenesum FITS file (default: ifs_im.scenesum.fits).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-o, --output%s %s<file>%s        "
           "Output ASCII file path (default: score.txt).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-h, --help%s                  "
           "Show this help message.\n\n",
           c_green, c_reset);

    printf("%sCOLOR MODE%s\n", c_cyan, c_reset);
    printf("  Colors: %s%s%s. Set %sNO_COLOR%s env var to disable.\n",
           c_yellow, (getenv("NO_COLOR") == NULL) ? "ENABLED" : "DISABLED", c_reset,
           c_bold, c_reset);
}

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

int compare_pixels(
    const void *a,
    const void *b)
{
    double val_a = ((PixelScore *)a)->val_planet;
    double val_b = ((PixelScore *)b)->val_planet;
    if (val_a > val_b) return -1;
    if (val_a < val_b) return 1;
    return 0;
}

int main(
    int   argc,
    char **argv)
{
    init_colors();

    char *planet_file   = "ifs_im.planet.fits";
    char *scenesum_file = "ifs_im.scenesum.fits";
    char *out_file      = "score.txt";

    static struct option long_options[] = {
        {"planet", required_argument, 0, 'p'},
        {"scenesum", required_argument, 0, 's'},
        {"output", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "p:s:o:h",
                             long_options, &option_index)) != -1)
    {
        switch (opt)
        {
            case 'p':
                planet_file = optarg;
                break;
            case 's':
                scenesum_file = optarg;
                break;
            case 'o':
                out_file = optarg;
                break;
            case 'h':
                print_help(argv[0]);
                return 0;
            default:
                print_help(argv[0]);
                return 1;
        }
    }

    /* Load pixel values */
    int p_rows = 0, p_cols = 0;
    double *planet_data = read_fits_2d(planet_file, &p_rows, &p_cols);
    if (planet_data == NULL)
    {
        fprintf(stderr, "ERROR: Failed to load planet FITS file: %s\n", planet_file);
        return 1;
    }

    int s_rows = 0, s_cols = 0;
    double *scenesum_data = read_fits_2d(scenesum_file, &s_rows, &s_cols);
    if (scenesum_data == NULL)
    {
        fprintf(stderr, "ERROR: Failed to load scenesum FITS file: %s\n", scenesum_file);
        free(planet_data);
        return 1;
    }

    if (p_rows != s_rows || p_cols != s_cols)
    {
        fprintf(stderr, "ERROR: Dimension mismatch. Planet is %dx%d, Scenesum is %dx%d.\n",
                p_rows, p_cols, s_rows, s_cols);
        free(planet_data);
        free(scenesum_data);
        return 1;
    }

    int M = p_rows;

    /* Extract WCS header coordinates from the planet FITS image */
    double cdelt1 = 0.5;
    double crpix1 = (M + 1) / 2.0;
    double crpix2 = (M + 1) / 2.0;

    fitsfile *fptr = NULL;
    int status = 0;
    fits_open_file(&fptr, planet_file, READONLY, &status);
    if (!status)
    {
        fits_read_key(fptr, TDOUBLE, "CDELT1", &cdelt1, NULL, &status);
        fits_read_key(fptr, TDOUBLE, "CRPIX1", &crpix1, NULL, &status);
        fits_read_key(fptr, TDOUBLE, "CRPIX2", &crpix2, NULL, &status);
        if (status)
        {
            fprintf(stderr, "WARNING: WCS keywords not fully present in %s. "
                            "Using standard grid fallback (CDELT=0.5).\n",
                            planet_file);
            cdelt1 = 0.5;
            crpix1 = (M + 1) / 2.0;
            crpix2 = (M + 1) / 2.0;
            status = 0;
        }
        fits_close_file(fptr, &status);
    }
    else
    {
        status = 0;
    }

    int total_pixels = M * M;
    PixelScore *pixels = malloc(total_pixels * sizeof(PixelScore));
    if (pixels == NULL)
    {
        fprintf(stderr, "ERROR: Out of memory for scoring calculations.\n");
        free(planet_data);
        free(scenesum_data);
        return 1;
    }

    for (int ly = 0; ly < M; ly++)
    {
        for (int lx = 0; lx < M; lx++)
        {
            int idx = ly * M + lx;
            PixelScore *p = &pixels[idx];
            p->row = ly;
            p->col = lx;
            p->x_ld = (lx + 1.0 - crpix1) * cdelt1;
            p->y_ld = (ly + 1.0 - crpix2) * cdelt1;
            p->val_planet = planet_data[idx];
            p->val_scenesum = scenesum_data[idx];
            p->val_diff = p->val_scenesum - p->val_planet;
        }
    }

    free(planet_data);
    free(scenesum_data);

    /* Sort by decreasing planet intensity */
    qsort(pixels, total_pixels, sizeof(PixelScore), compare_pixels);

    FILE *out = fopen(out_file, "w");
    if (out == NULL)
    {
        fprintf(stderr, "ERROR: Could not open output file %s\n", out_file);
        free(pixels);
        return 1;
    }

    /* Print Header */
    fprintf(out, "# MEFS Score Result File\n");
    fprintf(out, "# Planet FITS file:   %s\n", planet_file);
    fprintf(out, "# Scenesum FITS file: %s\n", scenesum_file);
    fprintf(out, "# Pixel scale (CDELT): %.4f lambda/D\n", cdelt1);
    fprintf(out, "# Ref Pixel (CRPIX):  (%.2f, %.2f)\n", crpix1, crpix2);
    fprintf(out, "# Format:\n");
    fprintf(out, "# Col 1:  Pixel Rank index (1-based, ordered by decreasing planet intensity)\n");
    fprintf(out, "# Col 2:  X coordinate (lambda/D)\n");
    fprintf(out, "# Col 3:  Y coordinate (lambda/D)\n");
    fprintf(out, "# Col 4:  Planet intensity (ifs_im.planet.fits)\n");
    fprintf(out, "# Col 5:  Total intensity (ifs_im.scenesum.fits)\n");
    fprintf(out, "# Col 6:  Difference (Total - Planet, non-planet light)\n");
    fprintf(out, "# Col 7:  Cumulative Planet intensity\n");
    fprintf(out, "# Col 8:  Cumulative Total intensity\n");
    fprintf(out, "# Col 9:  Cumulative Difference intensity\n");
    fprintf(out, "# Col 10: SNR for this pixel\n");
    fprintf(out, "# Col 11: Cumulative SNR\n");
    fprintf(out, "# Col 12: Row index\n");
    fprintf(out, "# Col 13: Col index\n");
    fprintf(out, "# ------------------------------------------------------------"
                 "--------------------\n");

    double cum_planet = 0.0;
    double cum_scenesum = 0.0;
    double cum_diff = 0.0;
    double cum_snr_sq = 0.0;

    for (int ii = 0; ii < total_pixels; ii++)
    {
        PixelScore *p = &pixels[ii];
        cum_planet += p->val_planet;
        cum_scenesum += p->val_scenesum;
        cum_diff += p->val_diff;

        double pixel_snr = 0.0;
        if (p->val_planet + p->val_scenesum > 0.0)
        {
            pixel_snr = p->val_planet / sqrt(p->val_planet + p->val_scenesum);
        }
        cum_snr_sq += pixel_snr * pixel_snr;
        double cum_snr = sqrt(cum_snr_sq);

        fprintf(out,
                "%7d %10.4f %10.4f %15.6e %15.6e %15.6e "
                "%15.6e %15.6e %15.6e %15.6e %15.6e %5d %5d\n",
                ii + 1, p->x_ld, p->y_ld, p->val_planet, p->val_scenesum, p->val_diff,
                cum_planet, cum_scenesum, cum_diff, pixel_snr, cum_snr, p->row, p->col);
    }

    fclose(out);
    free(pixels);
    printf("Successfully wrote score list to: %s\n", out_file);
    return 0;
}
