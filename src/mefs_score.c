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

enum
{
    OPT_PIXEL_SCALE = 256
};

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
    printf("  %s%s%s - Grade planet detection "
           "contrast metric from IFS outputs.\n\n",
           c_b_green, prog_name, c_reset);

    printf("%sUSAGE%s\n", c_cyan, c_reset);
    printf("  %s%s%s [OPTIONS]\n\n",
           c_b_green, prog_name, c_reset);

    printf("%sDESCRIPTION%s\n", c_cyan, c_reset);
    printf("  For each pixel in the planet image,"
           " measures total, planet, and\n");
    printf("  background flux, and outputs a "
           "sorted list by planet intensity.\n");
    printf("  The output ASCII file contains "
           "cumulative integrations to evaluate\n");
    printf("  detection metrics.\n\n");
    printf("  Two modes are supported:\n");
    printf("    %sCollapsed IFS%s (default): "
           "Uses lenslet-integrated flux maps\n",
           c_bold, c_reset);
    printf("    (%sifs_im.*.fits%s, M x M "
           "pixels). Pixel scale from WCS.\n",
           c_yellow, c_reset);
    printf("    %sRaw IFS%s (-r): Uses full-"
           "resolution detector images\n",
           c_bold, c_reset);
    printf("    (%sifsraw_im.*.fits%s, M*P x "
           "M*P pixels). Pixel scale\n",
           c_yellow, c_reset);
    printf("    computed from lenslet count and"
           " lenslet pixel scale.\n\n");

    printf("%sOPTIONS%s\n", c_cyan, c_reset);
    printf("  %s-p, --planet%s %s<file>%s"
           "        Planet flux map FITS file.\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("                       "
           "Default: ifs_im.planet.fits "
           "(or ifsraw_im.planet.fits with -r)\n");
    printf("  %s-s, --scenesum%s %s<file>%s"
           "      Combined scenesum FITS file.\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("                       "
           "Default: ifs_im.scenesum.fits "
           "(or ifsraw_im.scenesum.fits with -r)\n");
    printf("  %s-o, --output%s %s<file>%s"
           "        Output ASCII file path "
           "(default: score.txt).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("  %s-r, --raw%s               "
           "   Use raw IFS images instead of "
           "collapsed.\n",
           c_green, c_reset);
    printf("                       "
           "Default output becomes "
           "score_raw.txt.\n");
    printf("  %s-m, --lenslet_count%s %s<N>%s"
           "   Lenslet grid dimension "
           "(default: 32).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("                       "
           "Only used in raw mode to compute"
           " pixel scale.\n");
    printf("  %s--pixel_scale%s %s<val>%s"
           "       Lenslet pixel scale in "
           "lambda/D (default: 0.5).\n",
           c_green, c_reset, c_magenta, c_reset);
    printf("                       "
           "Only used in raw mode. "
           "Raw pixel scale = val / P.\n");
    printf("  %s-h, --help%s               "
           "   Show this help message.\n\n",
           c_green, c_reset);

    printf("%sEXAMPLES%s\n", c_cyan, c_reset);
    printf("  1. Score collapsed IFS "
           "(default):\n");
    printf("     %s$%s %s%s%s\n",
           c_grey, c_reset,
           c_b_green, prog_name, c_reset);
    printf("  2. Score raw IFS:\n");
    printf("     %s$%s %s%s -r%s\n",
           c_grey, c_reset,
           c_b_green, prog_name, c_reset);
    printf("  3. Score raw IFS with custom "
           "lenslet count:\n");
    printf("     %s$%s %s%s -r -m 64%s\n\n",
           c_grey, c_reset,
           c_b_green, prog_name, c_reset);

    printf("%sCOLOR MODE%s\n", c_cyan, c_reset);
    printf("  Colors: %s%s%s. Set %sNO_COLOR%s"
           " env var to disable.\n",
           c_yellow,
           (getenv("NO_COLOR") == NULL)
               ? "ENABLED" : "DISABLED",
           c_reset, c_bold, c_reset);
}

double *read_fits_2d(
    const char *filename,
    int        *out_rows,
    int        *out_cols)
{
    fitsfile *fptr = NULL;
    int status = 0;

    fits_open_file(&fptr, filename,
                   READONLY, &status);
    if (status)
    {
        fprintf(stderr,
                "ERROR: [%s:%d] Error opening "
                "FITS file %s: ",
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
        fprintf(stderr,
                "ERROR: [%s:%d] Error reading "
                "image dimensions for %s: ",
                __func__, __LINE__, filename);
        fits_report_error(stderr, status);
        fits_close_file(fptr, &status);
        return NULL;
    }

    *out_rows = naxes[1];
    *out_cols = naxes[0];

    long npix = naxes[0] * naxes[1];
    double *data = malloc(npix * sizeof(double));
    if (data == NULL)
    {
        fprintf(stderr,
                "ERROR: [%s:%d] Failed to "
                "allocate memory for image "
                "data.\n",
                __func__, __LINE__);
        fits_close_file(fptr, &status);
        return NULL;
    }

    long fpixel[2] = {1, 1};
    fits_read_pix(fptr, TDOUBLE, fpixel,
                  npix, NULL, data,
                  NULL, &status);
    if (status)
    {
        fprintf(stderr,
                "ERROR: [%s:%d] Error reading "
                "pixels from FITS: ",
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

    char *planet_file   = NULL;
    char *scenesum_file = NULL;
    char *out_file      = NULL;
    int raw_mode        = 0;
    int lenslet_count   = 32;
    double pixel_scale  = 0.5;

    /* Track whether user explicitly set files */
    int user_set_planet   = 0;
    int user_set_scenesum = 0;
    int user_set_output   = 0;

    static struct option long_options[] = {
        {"planet",        required_argument,
         0, 'p'},
        {"scenesum",      required_argument,
         0, 's'},
        {"output",        required_argument,
         0, 'o'},
        {"raw",           no_argument,
         0, 'r'},
        {"lenslet_count", required_argument,
         0, 'm'},
        {"pixel_scale",   required_argument,
         0, OPT_PIXEL_SCALE},
        {"help",          no_argument,
         0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(
                argc, argv, "p:s:o:rm:h",
                long_options,
                &option_index)) != -1)
    {
        switch (opt)
        {
            case 'p':
                planet_file = optarg;
                user_set_planet = 1;
                break;
            case 's':
                scenesum_file = optarg;
                user_set_scenesum = 1;
                break;
            case 'o':
                out_file = optarg;
                user_set_output = 1;
                break;
            case 'r':
                raw_mode = 1;
                break;
            case 'm':
                lenslet_count = atoi(optarg);
                if (lenslet_count <= 0)
                {
                    fprintf(stderr,
                            "ERROR: lenslet_count"
                            " must be positive.\n");
                    return 1;
                }
                break;
            case OPT_PIXEL_SCALE:
                pixel_scale = atof(optarg);
                if (pixel_scale <= 0.0)
                {
                    fprintf(stderr,
                            "ERROR: pixel_scale"
                            " must be positive.\n");
                    return 1;
                }
                break;
            case 'h':
                print_help(argv[0]);
                return 0;
            default:
                print_help(argv[0]);
                return 1;
        }
    }

    /* Apply mode-dependent defaults */
    if (raw_mode)
    {
        if (!user_set_planet)
        {
            planet_file = "ifsraw_im.planet.fits";
        }
        if (!user_set_scenesum)
        {
            scenesum_file =
                "ifsraw_im.scenesum.fits";
        }
        if (!user_set_output)
        {
            out_file = "score_raw.txt";
        }
        printf("Mode: %sRaw IFS%s "
               "(full-resolution detector "
               "pixels)\n",
               c_bold, c_reset);
    }
    else
    {
        if (!user_set_planet)
        {
            planet_file = "ifs_im.planet.fits";
        }
        if (!user_set_scenesum)
        {
            scenesum_file =
                "ifs_im.scenesum.fits";
        }
        if (!user_set_output)
        {
            out_file = "score.txt";
        }
        printf("Mode: %sCollapsed IFS%s "
               "(lenslet-integrated flux)\n",
               c_bold, c_reset);
    }

    /* Load pixel values */
    printf("Reading planet flux map:   "
           "%s%s%s\n",
           c_yellow, planet_file, c_reset);
    int p_rows = 0, p_cols = 0;
    double *planet_data = read_fits_2d(
        planet_file, &p_rows, &p_cols);
    if (planet_data == NULL)
    {
        fprintf(stderr,
                "ERROR: Failed to load planet "
                "FITS file: %s\n", planet_file);
        return 1;
    }
    printf("  -> %dx%d pixels\n",
           p_cols, p_rows);

    printf("Reading scenesum flux map: "
           "%s%s%s\n",
           c_yellow, scenesum_file, c_reset);
    int s_rows = 0, s_cols = 0;
    double *scenesum_data = read_fits_2d(
        scenesum_file, &s_rows, &s_cols);
    if (scenesum_data == NULL)
    {
        fprintf(stderr,
                "ERROR: Failed to load scenesum "
                "FITS file: %s\n",
                scenesum_file);
        free(planet_data);
        return 1;
    }
    printf("  -> %dx%d pixels\n",
           s_cols, s_rows);

    if (p_rows != s_rows || p_cols != s_cols)
    {
        fprintf(stderr,
                "ERROR: Dimension mismatch. "
                "Planet is %dx%d, "
                "Scenesum is %dx%d.\n",
                p_rows, p_cols,
                s_rows, s_cols);
        free(planet_data);
        free(scenesum_data);
        return 1;
    }

    int N = p_rows;

    /* Compute pixel scale and origin */
    double cdelt1 = 0.5;
    double crpix1 = (N + 1) / 2.0;
    double crpix2 = (N + 1) / 2.0;

    if (raw_mode)
    {
        /* Raw mode: derive pixel scale from
           lenslet count and lenslet pixel
           scale. P = image_size / M */
        if (N % lenslet_count != 0)
        {
            fprintf(stderr,
                    "ERROR: Image size %d is not "
                    "divisible by lenslet_count "
                    "%d.\n",
                    N, lenslet_count);
            free(planet_data);
            free(scenesum_data);
            return 1;
        }
        int P = N / lenslet_count;
        cdelt1 = pixel_scale / P;
        crpix1 = (N + 1) / 2.0;
        crpix2 = (N + 1) / 2.0;
        printf("Raw IFS pixel scale:      "
               "P=%d sub-pixels/lenslet\n", P);
        printf("  -> lenslet scale = %.4f l/D"
               "  raw pixel scale = %.6f l/D\n",
               pixel_scale, cdelt1);
        printf("  -> origin at (%.1f, %.1f)\n",
               crpix1, crpix2);
    }
    else
    {
        /* Collapsed mode: extract WCS from
           planet FITS header */
        printf("Extracting WCS from:      "
               "%s%s%s\n",
               c_yellow, planet_file, c_reset);

        fitsfile *fptr = NULL;
        int status = 0;
        fits_open_file(&fptr, planet_file,
                       READONLY, &status);
        if (!status)
        {
            fits_read_key(fptr, TDOUBLE,
                          "CDELT1", &cdelt1,
                          NULL, &status);
            fits_read_key(fptr, TDOUBLE,
                          "CRPIX1", &crpix1,
                          NULL, &status);
            fits_read_key(fptr, TDOUBLE,
                          "CRPIX2", &crpix2,
                          NULL, &status);
            if (status)
            {
                fprintf(stderr,
                        "WARNING: WCS keywords "
                        "not fully present in "
                        "%s. Using standard grid"
                        " fallback "
                        "(CDELT=0.5).\n",
                        planet_file);
                cdelt1 = 0.5;
                crpix1 = (N + 1) / 2.0;
                crpix2 = (N + 1) / 2.0;
                status = 0;
            }
            fits_close_file(fptr, &status);
            printf("  -> CDELT1=%.4f l/D  "
                   "CRPIX=(%.1f, %.1f)\n",
                   cdelt1, crpix1, crpix2);
        }
        else
        {
            status = 0;
        }
    } // if (raw_mode)

    int total_pixels = N * N;
    PixelScore *pixels =
        malloc(total_pixels * sizeof(PixelScore));
    if (pixels == NULL)
    {
        fprintf(stderr,
                "ERROR: Out of memory for "
                "scoring calculations.\n");
        free(planet_data);
        free(scenesum_data);
        return 1;
    }

    for (int ly = 0; ly < N; ly++)
    {
        for (int lx = 0; lx < N; lx++)
        {
            int idx = ly * N + lx;
            PixelScore *p = &pixels[idx];
            p->row = ly;
            p->col = lx;
            p->x_ld =
                (lx + 1.0 - crpix1) * cdelt1;
            p->y_ld =
                (ly + 1.0 - crpix2) * cdelt1;
            p->val_planet = planet_data[idx];
            p->val_scenesum =
                scenesum_data[idx];
            p->val_diff =
                p->val_scenesum - p->val_planet;
        }
    }

    free(planet_data);
    free(scenesum_data);

    /* Sort by decreasing planet intensity */
    printf("Sorting %d pixels by planet "
           "intensity...\n", total_pixels);
    qsort(pixels, total_pixels,
          sizeof(PixelScore), compare_pixels);

    printf("Writing score output:     "
           "%s%s%s\n",
           c_yellow, out_file, c_reset);
    FILE *out = fopen(out_file, "w");
    if (out == NULL)
    {
        fprintf(stderr,
                "ERROR: Could not open output "
                "file %s\n", out_file);
        free(pixels);
        return 1;
    }

    /* Print Header */
    fprintf(out,
            "# MEFS Score Result File\n");
    fprintf(out,
            "# Mode: %s\n",
            raw_mode ? "Raw IFS" : "Collapsed IFS");
    fprintf(out,
            "# Planet FITS file:   %s\n",
            planet_file);
    fprintf(out,
            "# Scenesum FITS file: %s\n",
            scenesum_file);
    fprintf(out,
            "# Pixel scale (CDELT): "
            "%.6f lambda/D\n", cdelt1);
    fprintf(out,
            "# Ref Pixel (CRPIX):  "
            "(%.2f, %.2f)\n", crpix1, crpix2);
    if (raw_mode)
    {
        int P = N / lenslet_count;
        fprintf(out,
                "# Lenslet count: %d  "
                "Sub-pixels/lenslet: %d\n",
                lenslet_count, P);
        fprintf(out,
                "# Lenslet pixel scale: "
                "%.4f lambda/D\n", pixel_scale);
    }
    fprintf(out, "# Format:\n");
    fprintf(out,
            "# Col 1:  Pixel Rank index "
            "(1-based, ordered by decreasing "
            "planet intensity)\n");
    fprintf(out,
            "# Col 2:  X coordinate "
            "(lambda/D)\n");
    fprintf(out,
            "# Col 3:  Y coordinate "
            "(lambda/D)\n");
    fprintf(out,
            "# Col 4:  Planet intensity "
            "(%s)\n", planet_file);
    fprintf(out,
            "# Col 5:  Total intensity "
            "(%s)\n", scenesum_file);
    fprintf(out,
            "# Col 6:  Difference "
            "(Total - Planet, non-planet "
            "light)\n");
    fprintf(out,
            "# Col 7:  Cumulative Planet "
            "intensity\n");
    fprintf(out,
            "# Col 8:  Cumulative Total "
            "intensity\n");
    fprintf(out,
            "# Col 9:  Cumulative Difference "
            "intensity\n");
    fprintf(out,
            "# Col 10: SNR for this pixel\n");
    fprintf(out,
            "# Col 11: Cumulative SNR\n");
    fprintf(out,
            "# Col 12: Row index\n");
    fprintf(out,
            "# Col 13: Col index\n");
    fprintf(out,
            "# ---------------------------"
            "----------------------------"
            "-------------------------\n");

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
        if (p->val_planet + p->val_scenesum
            > 0.0)
        {
            pixel_snr = p->val_planet
                / sqrt(p->val_planet
                       + p->val_scenesum);
        }
        cum_snr_sq += pixel_snr * pixel_snr;
        double cum_snr = sqrt(cum_snr_sq);

        fprintf(out,
                "%7d %10.4f %10.4f "
                "%15.6e %15.6e %15.6e "
                "%15.6e %15.6e %15.6e "
                "%15.6e %15.6e %5d %5d\n",
                ii + 1,
                p->x_ld, p->y_ld,
                p->val_planet,
                p->val_scenesum,
                p->val_diff,
                cum_planet,
                cum_scenesum,
                cum_diff,
                pixel_snr, cum_snr,
                p->row, p->col);
    }

    fclose(out);
    free(pixels);
    printf("Successfully wrote score list "
           "to: %s\n", out_file);
    return 0;
}
