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
    OPT_PIXEL_SCALE = 256,
    OPT_DARK        = 257,
    OPT_FLUX_SCALE  = 258,
    OPT_SWEEP       = 259
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
    printf("  For each pixel in the planet "
           "image, measures total, planet, and\n");
    printf("  background flux, and outputs a "
           "sorted list by planet intensity.\n");
    printf("  The output ASCII file contains "
           "cumulative integrations to evaluate\n");
    printf("  detection metrics. SNR values "
           "assume a 1-second exposure;\n");
    printf("  SNR scales as sqrt(t). "
           "Exposure time to reach\n");
    printf("  SNR=20 is reported in the "
           "summary.\n\n");
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
    printf("    computed from lenslet count"
           " and lenslet pixel scale.\n\n");

    printf("  If %smefs.log%s is present in "
           "the current directory, the program\n",
           c_yellow, c_reset);
    printf("  also computes MEFS SNR (with "
           "dark current) and compares IFS\n");
    printf("  vs MEFS detection performance."
           "\n\n");

    printf("%sOPTIONS%s\n", c_cyan, c_reset);
    printf("  %s-p, --planet%s %s<file>%s"
           "        Planet flux map FITS "
           "file.\n",
           c_green, c_reset,
           c_magenta, c_reset);
    printf("                       "
           "Default: ifs_im.planet.fits "
           "(or ifsraw_im.planet.fits "
           "with -r)\n");
    printf("  %s-s, --scenesum%s %s<file>%s"
           "      Combined scenesum FITS "
           "file.\n",
           c_green, c_reset,
           c_magenta, c_reset);
    printf("                       "
           "Default: ifs_im.scenesum.fits "
           "(or ifsraw_im.scenesum.fits "
           "with -r)\n");
    printf("  %s-o, --output%s %s<file>%s"
           "        Output ASCII file path "
           "(default: score.txt).\n",
           c_green, c_reset,
           c_magenta, c_reset);
    printf("  %s-r, --raw%s               "
           "   Use raw IFS images instead"
           " of collapsed.\n",
           c_green, c_reset);
    printf("                       "
           "Default output becomes "
           "score_raw.txt.\n");
    printf("  %s-d, --dark%s %s<val>%s"
           "          Dark current per "
           "pixel (default: 0.0).\n",
           c_green, c_reset,
           c_magenta, c_reset);
    printf("                       "
           "Added as incoherent noise "
           "under sqrt() in SNR.\n");
    printf("  %s-m, --lenslet_count%s "
           "%s<N>%s   Lenslet grid "
           "dimension (default: 32).\n",
           c_green, c_reset,
           c_magenta, c_reset);
    printf("                       "
           "Only used in raw mode to "
           "compute pixel scale.\n");
    printf("  %s-F, --flux_scale%s %s<val>%s"
           "    Flux scaling factor for"
           " all sources (default: "
           "1.0).\n",
           c_green, c_reset,
           c_magenta, c_reset);
    printf("                       "
           "Scales planet, star, zodi,"
           " exozodi. Dark current is"
           " not scaled.\n");
    printf("  %s--pixel_scale%s %s<val>%s"
           "       Lenslet pixel scale"
           " in lambda/D (default: "
           "0.5).\n",
           c_green, c_reset,
           c_magenta, c_reset);
    printf("                       "
           "Only used in raw mode. "
           "Raw pixel scale = val / P."
           "\n");
    printf("  %s--sweep%s                 "
           "   Sweep flux multiplier R"
           " from 10 to 1e5\n",
           c_green, c_reset);
    printf("                       "
           "(steps of 1.1x). Writes "
           "%sscore_sweep.txt%s with\n",
           c_yellow, c_reset);
    printf("                       "
           "IFS and MEFS exposure "
           "times for SNR=20.\n");
    printf("  %s-h, --help%s               "
           "   Show this help message."
           "\n\n",
           c_green, c_reset);

    printf("%sEXAMPLES%s\n", c_cyan, c_reset);
    printf("  1. Score collapsed IFS "
           "(default):\n");
    printf("     %s$%s %s%s%s\n",
           c_grey, c_reset,
           c_b_green, prog_name, c_reset);
    printf("  2. Score with dark current:\n");
    printf("     %s$%s %s%s --dark 0.001%s\n",
           c_grey, c_reset,
           c_b_green, prog_name, c_reset);
    printf("  3. Score raw IFS:\n");
    printf("     %s$%s %s%s -r%s\n",
           c_grey, c_reset,
           c_b_green, prog_name, c_reset);
    printf("  4. Score raw IFS with custom "
           "lenslet count:\n");
    printf("     %s$%s %s%s -r -m 64%s\n\n",
           c_grey, c_reset,
           c_b_green, prog_name, c_reset);

    printf("%sMEFS COMPARISON%s\n",
           c_cyan, c_reset);
    printf("  If %smefs.log%s exists, the "
           "program reads MEFS coupling\n",
           c_yellow, c_reset);
    printf("  data and computes MEFS SNR "
           "with the specified dark\n");
    printf("  current (single pixel). A "
           "summary comparing IFS and\n");
    printf("  MEFS SNR is appended to the "
           "output file.\n\n");

    printf("%sCOLOR MODE%s\n",
           c_cyan, c_reset);
    printf("  Colors: %s%s%s. Set %sNO_COLOR"
           "%s env var to disable.\n",
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

/**
 * parse_mefs_log - Read MEFS coupling data
 *                  from mefs.log
 * @mefs_planet:    [out] planet coupled light
 * @mefs_background:[out] background coupled light
 * @mefs_snr_orig:  [out] original MEFS SNR
 *
 * Return: 0 on success, -1 if file not found or
 *         parse error.
 */
static int parse_mefs_log(
    double *mefs_planet,
    double *mefs_background,
    double *mefs_snr_orig)
{
    FILE *f = fopen("mefs.log", "r");
    if (f == NULL)
    {
        return -1;
    }

    *mefs_planet     = 0.0;
    *mefs_background = 0.0;
    *mefs_snr_orig   = 0.0;

    int found_planet = 0;
    int found_bg     = 0;
    int found_snr    = 0;

    char line[512];
    while (fgets(line, sizeof(line), f) != NULL)
    {
        double val;
        if (sscanf(line,
                   "Planet Light Gathered: %le",
                   &val) == 1)
        {
            *mefs_planet = val;
            found_planet = 1;
        }
        else if (sscanf(line,
                        "Total Background "
                        "Light: %le",
                        &val) == 1)
        {
            *mefs_background = val;
            found_bg = 1;
        }
        else if (sscanf(line,
                        "Combined Spectrograph "
                        "SNR: %le",
                        &val) == 1)
        {
            *mefs_snr_orig = val;
            found_snr = 1;
        }
    } // while reading lines

    fclose(f);

    if (!found_planet || !found_bg || !found_snr)
    {
        return -1;
    }
    return 0;
}

int compare_pixels(
    const void *a,
    const void *b)
{
    double val_a =
        ((PixelScore *)a)->val_planet;
    double val_b =
        ((PixelScore *)b)->val_planet;
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
    double dark_current = 0.0;
    double flux_scale   = 1.0;
    int sweep_mode      = 0;

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
        {"dark",          required_argument,
         0, 'd'},
        {"lenslet_count", required_argument,
         0, 'm'},
        {"pixel_scale",   required_argument,
         0, OPT_PIXEL_SCALE},
        {"flux_scale",    required_argument,
         0, OPT_FLUX_SCALE},
        {"sweep",         no_argument,
         0, OPT_SWEEP},
        {"help",          no_argument,
         0, 'h'},
        {0, 0, 0, 0}
    };

    int opt;
    int option_index = 0;
    while ((opt = getopt_long(
                argc, argv, "p:s:o:rd:m:F:h",
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
            case 'd':
                dark_current = atof(optarg);
                if (dark_current < 0.0)
                {
                    fprintf(stderr,
                            "ERROR: dark current"
                            " must be >= 0.\n");
                    return 1;
                }
                break;
            case 'm':
                lenslet_count = atoi(optarg);
                if (lenslet_count <= 0)
                {
                    fprintf(stderr,
                            "ERROR: lenslet_count"
                            " must be positive."
                            "\n");
                    return 1;
                }
                break;
            case OPT_PIXEL_SCALE:
                pixel_scale = atof(optarg);
                if (pixel_scale <= 0.0)
                {
                    fprintf(stderr,
                            "ERROR: pixel_scale"
                            " must be positive."
                            "\n");
                    return 1;
                }
                break;
            case OPT_FLUX_SCALE:
            case 'F':
                flux_scale = atof(optarg);
                if (flux_scale <= 0.0)
                {
                    fprintf(stderr,
                            "ERROR: flux_scale"
                            " must be "
                            "positive.\n");
                    return 1;
                }
                break;
            case OPT_SWEEP:
                sweep_mode = 1;
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
            planet_file =
                "ifsraw_im.planet.fits";
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
            planet_file =
                "ifs_im.planet.fits";
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

    if (dark_current > 0.0)
    {
        printf("Dark current per pixel: "
               "%s%.6e%s\n",
               c_yellow, dark_current,
               c_reset);
    }
    if (flux_scale != 1.0)
    {
        printf("Flux scaling factor:    "
               "%s%.6e%s\n",
               c_yellow, flux_scale,
               c_reset);
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
                    "ERROR: Image size %d is "
                    "not divisible by "
                    "lenslet_count %d.\n",
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
        printf("  -> lenslet scale = %.4f "
               "l/D  raw pixel scale = "
               "%.6f l/D\n",
               pixel_scale, cdelt1);
        printf("  -> origin at (%.1f, %.1f)"
               "\n", crpix1, crpix2);
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
                        "%s. Using standard "
                        "grid fallback "
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
        malloc(total_pixels
               * sizeof(PixelScore));
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
                p->val_scenesum
                - p->val_planet;
        }
    }

    /* Apply flux scaling factor */
    if (flux_scale != 1.0)
    {
        for (int ii = 0;
             ii < total_pixels; ii++)
        {
            pixels[ii].val_planet *=
                flux_scale;
            pixels[ii].val_scenesum *=
                flux_scale;
            pixels[ii].val_diff *=
                flux_scale;
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
            raw_mode
                ? "Raw IFS"
                : "Collapsed IFS");
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
            "(%.2f, %.2f)\n",
            crpix1, crpix2);
    fprintf(out,
            "# Dark current/pixel: "
            "%.6e\n", dark_current);
    fprintf(out,
            "# Flux scaling factor: "
            "%.6e\n", flux_scale);
    if (raw_mode)
    {
        int P = N / lenslet_count;
        fprintf(out,
                "# Lenslet count: %d  "
                "Sub-pixels/lenslet: %d\n",
                lenslet_count, P);
        fprintf(out,
                "# Lenslet pixel scale: "
                "%.4f lambda/D\n",
                pixel_scale);
    }
    fprintf(out, "# Format:\n");
    fprintf(out,
            "# Col 1:  Pixel Rank index "
            "(1-based, decreasing planet "
            "intensity)\n");
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
            "(Total - Planet)\n");
    fprintf(out,
            "# Col 7:  Cumulative Planet "
            "intensity\n");
    fprintf(out,
            "# Col 8:  Cumulative Total "
            "intensity\n");
    fprintf(out,
            "# Col 9:  Cumulative "
            "Difference intensity\n");
    fprintf(out,
            "# Col 10: Cumulative dark "
            "current (rank * dark)\n");
    fprintf(out,
            "# Col 11: SNR for this "
            "pixel\n");
    fprintf(out,
            "# Col 12: Cumulative SNR\n");
    fprintf(out,
            "# Col 13: Row index\n");
    fprintf(out,
            "# Col 14: Col index\n");
    fprintf(out,
            "# -------------------------"
            "---------------------------"
            "---------------------------"
            "\n");

    double cum_planet   = 0.0;
    double cum_scenesum = 0.0;
    double cum_diff     = 0.0;
    double cum_snr_sq   = 0.0;

    /* Track best cumulative SNR and its
       pixel count for the summary */
    double best_cum_snr = 0.0;
    int best_cum_npix   = 0;

    for (int ii = 0; ii < total_pixels; ii++)
    {
        PixelScore *p = &pixels[ii];
        cum_planet += p->val_planet;
        cum_scenesum += p->val_scenesum;
        cum_diff += p->val_diff;

        /* Per-pixel SNR with dark current:
           SNR = S_p / sqrt(S_p + S_total + D)
        */
        double denom =
            p->val_planet
            + p->val_scenesum
            + dark_current;
        double pixel_snr = 0.0;
        if (denom > 0.0)
        {
            pixel_snr =
                p->val_planet / sqrt(denom);
        }
        cum_snr_sq +=
            pixel_snr * pixel_snr;
        double cum_snr = sqrt(cum_snr_sq);

        /* Cumulative dark current */
        double cum_dark =
            (ii + 1) * dark_current;

        if (cum_snr > best_cum_snr)
        {
            best_cum_snr = cum_snr;
            best_cum_npix = ii + 1;
        }

        fprintf(out,
                "%7d %10.4f %10.4f "
                "%15.6e %15.6e %15.6e "
                "%15.6e %15.6e %15.6e "
                "%15.6e %15.6e %15.6e "
                "%5d %5d\n",
                ii + 1,
                p->x_ld, p->y_ld,
                p->val_planet,
                p->val_scenesum,
                p->val_diff,
                cum_planet,
                cum_scenesum,
                cum_diff,
                cum_dark,
                pixel_snr, cum_snr,
                p->row, p->col);
    } // for each pixel

    /* --- MEFS comparison summary --- */
    double mefs_planet   = 0.0;
    double mefs_bg       = 0.0;
    double mefs_snr_orig = 0.0;
    int have_mefs = (parse_mefs_log(
        &mefs_planet, &mefs_bg,
        &mefs_snr_orig) == 0);

    fprintf(out,
            "# \n"
            "# --- SNR COMPARISON ---\n");
    fprintf(out,
            "# Flux scaling factor:    "
            "%.6e\n", flux_scale);
    fprintf(out,
            "# Dark current per pixel: "
            "%.6e\n", dark_current);
    fprintf(out,
            "# IFS optimal cumulative SNR: "
            "%.6e  (N=%d pixels)\n",
            best_cum_snr, best_cum_npix);

    /* Exposure time to reach SNR=20:
       SNR(t) = SNR_1s * sqrt(t)
       t = (target / SNR_1s)^2 */
    double snr_target = 20.0;
    double t_ifs = 0.0;
    if (best_cum_snr > 0.0)
    {
        double r = snr_target / best_cum_snr;
        t_ifs = r * r;
    }
    fprintf(out,
            "# IFS time for SNR=%.0f:       "
            "%.2f sec (%.2f hr)\n",
            snr_target, t_ifs,
            t_ifs / 3600.0);

    /* Print to stdout as well */
    printf("\n%s--- SNR COMPARISON ---%s\n",
           c_cyan, c_reset);
    printf("Flux scaling factor:      "
           "%.6e\n", flux_scale);
    printf("Dark current per pixel:   "
           "%.6e\n", dark_current);
    printf("IFS optimal cumul. SNR:   "
           "%s%.6e%s  (N=%d pixels)\n",
           c_b_green, best_cum_snr,
           c_reset, best_cum_npix);
    printf("IFS time for SNR=%.0f:     "
           "%s%.2f sec (%.2f hr)%s\n",
           snr_target,
           c_yellow, t_ifs,
           t_ifs / 3600.0, c_reset);

    if (have_mefs)
    {
        printf("Reading MEFS data from:   "
               "%smefs.log%s\n",
               c_yellow, c_reset);

        /* Scale MEFS values by flux factor */
        mefs_planet *= flux_scale;
        mefs_bg     *= flux_scale;

        /* MEFS SNR with dark current:
           all coupled light on 1 pixel */
        double mefs_denom =
            mefs_planet + mefs_bg
            + dark_current;
        double mefs_snr_dark = 0.0;
        if (mefs_denom > 0.0)
        {
            mefs_snr_dark =
                mefs_planet
                / sqrt(mefs_denom);
        }

        double ratio = 0.0;
        if (best_cum_snr > 0.0)
        {
            ratio =
                mefs_snr_dark / best_cum_snr;
        }

        fprintf(out,
                "# MEFS planet coupled:    "
                "%.6e\n", mefs_planet);
        fprintf(out,
                "# MEFS background coupled:"
                " %.6e\n", mefs_bg);
        fprintf(out,
                "# MEFS SNR (no dark):     "
                "%.6e\n", mefs_snr_orig);
        fprintf(out,
                "# MEFS SNR (with dark):   "
                "%.6e\n", mefs_snr_dark);
        fprintf(out,
                "# MEFS / IFS SNR ratio:   "
                "%.6f\n", ratio);

        double t_mefs = 0.0;
        if (mefs_snr_dark > 0.0)
        {
            double r =
                snr_target / mefs_snr_dark;
            t_mefs = r * r;
        }
        fprintf(out,
                "# MEFS time for SNR=%.0f:  "
                "    %.2f sec (%.2f hr)\n",
                snr_target, t_mefs,
                t_mefs / 3600.0);
        fprintf(out,
                "# IFS/MEFS time ratio:    "
                "%.6f\n",
                (t_mefs > 0.0)
                    ? t_ifs / t_mefs : 0.0);

        printf("  MEFS planet coupled:    "
               "%.6e\n", mefs_planet);
        printf("  MEFS background coupled:"
               " %.6e\n", mefs_bg);
        printf("  MEFS SNR (no dark):     "
               "%.6e\n", mefs_snr_orig);
        printf("  MEFS SNR (with dark):   "
               "%s%.6e%s\n",
               c_b_green, mefs_snr_dark,
               c_reset);
        printf("  MEFS / IFS SNR ratio:   "
               "%s%.6f%s\n",
               c_yellow, ratio, c_reset);

        printf("  MEFS time for SNR=%.0f:  "
               "  %s%.2f sec (%.2f hr)%s\n",
               snr_target,
               c_yellow, t_mefs,
               t_mefs / 3600.0, c_reset);
        printf("  IFS/MEFS time ratio:    "
               "%s%.6f%s\n",
               c_yellow,
               (t_mefs > 0.0)
                   ? t_ifs / t_mefs : 0.0,
               c_reset);
    }
    else
    {
        fprintf(out,
                "# MEFS: mefs.log not found"
                " (no comparison)\n");
        printf("MEFS: %smefs.log not found%s"
               " (no comparison)\n",
               c_grey, c_reset);
    }

    fclose(out);
    printf("\nSuccessfully wrote score list "
           "to: %s\n", out_file);

    /* --- Flux sweep mode --- */
    if (sweep_mode)
    {
        const char *sweep_file =
            "score_sweep.txt";
        FILE *sf = fopen(sweep_file, "w");
        if (sf == NULL)
        {
            fprintf(stderr,
                    "ERROR: Could not open "
                    "sweep file %s\n",
                    sweep_file);
            free(pixels);
            return 1;
        }

        fprintf(sf,
                "# MEFS Flux Sweep\n");
        fprintf(sf,
                "# Base flux_scale F = "
                "%.6e\n", flux_scale);
        fprintf(sf,
                "# Dark current/pixel = "
                "%.6e\n", dark_current);
        fprintf(sf,
                "# R sweeps from 10 to "
                "1e5, steps of 1.1x\n");
        fprintf(sf,
                "# Flux at each step = "
                "F / R  (pixel values "
                "already include F)\n");
        fprintf(sf,
                "# Col 1: R (divisor "
                "applied to F)\n");
        fprintf(sf,
                "# Col 2: Total flux "
                "scale (F/R)\n");
        fprintf(sf,
                "# Col 3: IFS optimal N "
                "(pixels)\n");
        fprintf(sf,
                "# Col 4: IFS SNR "
                "(1-sec)\n");
        fprintf(sf,
                "# Col 5: IFS time for "
                "SNR=20 (sec)\n");
        fprintf(sf,
                "# Col 6: IFS time for "
                "SNR=20 (hr)\n");
        if (have_mefs)
        {
            fprintf(sf,
                    "# Col 7: MEFS SNR "
                    "(1-sec)\n");
            fprintf(sf,
                    "# Col 8: MEFS time "
                    "for SNR=20 (sec)\n");
            fprintf(sf,
                    "# Col 9: MEFS time "
                    "for SNR=20 (hr)\n");
            fprintf(sf,
                    "# Col 10: MEFS/IFS "
                    "SNR ratio\n");
        }
        fprintf(sf,
                "# ----------------------"
                "----------------------"
                "----------------------"
                "\n");

        printf("\n%s--- FLUX SWEEP ---%s\n",
               c_cyan, c_reset);
        printf("Sweeping R from 10 to 1e5 "
               "(1.1x steps)...\n");

        /* Pixel values already have
           flux_scale applied. Sweep
           divides by R on top of that. */
        for (double R = 10.0;
             R <= 1.0e5 + 1.0;
             R *= 1.1)
        {
            /* IFS: find optimal cumulative
               SNR at this R */
            double best_snr = 0.0;
            int best_n = 0;
            double csq = 0.0;

            for (int ii = 0;
                 ii < total_pixels; ii++)
            {
                PixelScore *p = &pixels[ii];
                double sp = p->val_planet / R;
                double st =
                    p->val_scenesum / R;
                double den =
                    sp + st + dark_current;
                if (den > 0.0)
                {
                    double sn = sp / sqrt(den);
                    csq += sn * sn;
                }
                double cs = sqrt(csq);
                if (cs > best_snr)
                {
                    best_snr = cs;
                    best_n = ii + 1;
                }
            } // for each pixel

            /* IFS exposure time */
            double t_i = 0.0;
            if (best_snr > 0.0)
            {
                double rr =
                    snr_target / best_snr;
                t_i = rr * rr;
            }

            fprintf(sf,
                    "%15.6e %15.6e "
                    "%7d %15.6e "
                    "%15.2f %12.4f",
                    R, flux_scale / R,
                    best_n, best_snr,
                    t_i, t_i / 3600.0);

            if (have_mefs)
            {
                /* MEFS: single pixel,
                   scale coupled values
                   by 1/R (already scaled
                   by flux_scale above) */
                double mp = mefs_planet / R;
                double mb = mefs_bg / R;
                double md =
                    mp + mb + dark_current;
                double ms = 0.0;
                if (md > 0.0)
                {
                    ms = mp / sqrt(md);
                }
                double t_m = 0.0;
                if (ms > 0.0)
                {
                    double rr =
                        snr_target / ms;
                    t_m = rr * rr;
                }
                double ratio_s = 0.0;
                if (best_snr > 0.0)
                {
                    ratio_s = ms / best_snr;
                }
                fprintf(sf,
                        " %15.6e "
                        "%15.2f %12.4f "
                        "%12.6f",
                        ms, t_m,
                        t_m / 3600.0,
                        ratio_s);
            }
            fprintf(sf, "\n");
        } // for each R

        fclose(sf);
        printf("Wrote sweep to: "
               "%s%s%s\n",
               c_yellow, sweep_file,
               c_reset);
    } // if (sweep_mode)

    free(pixels);
    return 0;
}
