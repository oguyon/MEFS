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



#include "mefs_utils.h"
#include "mefs_fits.h"
#include "mefs_optics.h"
#include "mefs_simulation.h"

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
