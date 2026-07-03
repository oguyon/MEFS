#include "mefs_simulation.h"
#include "mefs_utils.h"
#include "mefs_fits.h"
#include "mefs_optics.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <glob.h>
#include <fftw3.h>
#ifdef _OPENMP
#include <omp.h>
#endif

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

